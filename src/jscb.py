#!/usr/bin/env python3
"""
One-to-one Python implementation of the KwanLab/JSCB pipeline:
- Reimplements the core Fortran clustering algorithm from jscb.f
- Reimplements the GI-calling logic from jscb.py

Input: GenBank file
Output: JSCB_output.gi (per gene cluster assignment) and JSCB_output.tsv (GI intervals)

Notes:
- This follows the repository behavior:
  * non-ATCG bases in genome are replaced with a random base
  * only CDS whose (end-start+1) % 3 == 0 are used
  * genes shorter than (mingenelen+3) codons are filtered out (matches Fortran logic)
  * "native cluster" = largest cluster by number of genes
  * GI = contiguous run of non-native-cluster genes of length >= 8 genes
  * reject GI if it fully contains an rRNA feature (same test as jscb.py)
"""

from __future__ import annotations

import math
import random
import sys
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
import os

from Bio import SeqIO


# -------------------------
# Math utilities (gamma)
# -------------------------

def _gammln(xx: float) -> float:
	# Port of Numerical Recipes gammln as in Fortran
	cof = [
		76.18009172947146,
		-86.50532032941677,
		24.01409824083091,
		-1.231739572450155,
		0.1208650973866179e-2,
		-0.5395239384953e-5,
	]
	stp = 2.5066282746310005
	x = xx
	y = x
	tmp = x + 5.5
	tmp = (x + 0.5) * math.log(tmp) - tmp
	ser = 1.000000000190015
	for c in cof:
		y += 1.0
		ser += c / y
	return tmp + math.log(stp * ser / x)


def _gser(a: float, x: float, itmax: int = 9000, eps: float = 3.0e-7) -> float:
	# Series representation of P(a, x)
	gln = _gammln(a)
	if x <= 0.0:
		return 0.0
	ap = a
	summ = 1.0 / a
	delt = summ
	for _ in range(1, itmax + 1):
		ap += 1.0
		delt *= x / ap
		summ += delt
		if abs(delt) < abs(summ) * eps:
			break
	return summ * math.exp(-x + a * math.log(x) - gln)


def _gcf(a: float, x: float, itmax: int = 9000, eps: float = 3.0e-7, fpmin: float = 1.0e-30) -> float:
	# Continued fraction representation of Q(a, x) = 1 - P(a, x)
	gln = _gammln(a)
	b = x + 1.0 - a
	c = 1.0 / fpmin
	d = 1.0 / b
	h = d
	for i in range(1, itmax + 1):
		an = -i * (i - a)
		b += 2.0
		d = an * d + b
		if abs(d) < fpmin:
			d = fpmin
		c = b + an / c
		if abs(c) < fpmin:
			c = fpmin
		d = 1.0 / d
		delt = d * c
		h *= delt
		if abs(delt - 1.0) < eps:
			break
	return math.exp(-x + a * math.log(x) - gln) * h


def gammp(a: float, x: float) -> float:
	"""
	Regularized lower incomplete gamma P(a, x), matching the Fortran gammp.
	"""
	if x < 0.0 or a <= 0.0:
		raise ValueError("bad arguments in gammp")
	if x < a + 1.0:
		return _gser(a, x)
	else:
		# gammp = 1 - gammcf
		return 1.0 - _gcf(a, x)


def sig(rnn: float, djsmx: float) -> float:
	"""
	Fortran:
	  snn = rnn*log(2)*djsmx
	  sig = gammp(30, snn)
	"""
	snn = rnn * math.log(2.0) * djsmx
	return gammp(30.0, snn)


# -------------------------
# Genetic code mapping from jscb.f (iaa / iab)
# -------------------------

def build_iaa_iab() -> Tuple[List[List[int]], List[int]]:
	"""
	Returns:
	  iaa: 20 x 6 (1-based codonIndex values), 0 where unused
	  iab: length 20, count of synonymous codons in that amino acid group
	"""
	iaa = [[0] * 6 for _ in range(20)]
	# These are direct ports of the Fortran assignments (lines 52..112)
	def set_row(r: int, vals: List[int]) -> None:
		for c, v in enumerate(vals):
			iaa[r - 1][c] = v

	set_row(1,  [42, 43])
	set_row(2,  [41, 44, 58, 59, 57, 60])
	set_row(3,  [26, 27, 25])
	set_row(4,  [28])
	set_row(5,  [74, 75, 73, 76])
	set_row(6,  [46, 47, 45, 48, 34, 35])
	set_row(7,  [62, 63, 61, 64])
	set_row(8,  [30, 31, 29, 32])
	set_row(9,  [78, 79, 77, 80])
	set_row(10, [38, 39])
	set_row(11, [54, 55])
	set_row(12, [53, 56])
	set_row(13, [22, 23])
	set_row(14, [21, 24])
	set_row(15, [70, 71])
	set_row(16, [69, 72])
	set_row(17, [50, 51])
	set_row(18, [52])
	set_row(19, [66, 67, 65, 68, 33, 36])
	set_row(20, [82, 83, 81, 84])

	# iab assignments (lines 114..133)
	iab = [
		2, 6, 3, 1, 4, 6, 4, 4, 4, 2,
		2, 2, 2, 2, 2, 2, 2, 1, 6, 4
	]
	return iaa, iab


# -------------------------
# Data structures
# -------------------------

@dataclass
class Gene:
	strand: str  # '+' or '-'
	start: int   # 1-based inclusive
	end: int	 # 1-based inclusive

	@property
	def length_nt(self) -> int:
		return self.end - self.start + 1


# -------------------------
# Codon indexing (matches Fortran logic)
# -------------------------

_BASE_TO_IY = {'A': 1, 'T': 2, 'C': 3, 'G': 4}
# complement mapping in terms of iy values used in Fortran '-' strand logic:
# if iy==1(A)->2(T), 2(T)->1(A), 3(C)->4(G), 4(G)->3(C)
_COMP_IY = {1: 2, 2: 1, 3: 4, 4: 3}


def genome_to_iy(genome: str) -> List[int]:
	# Fortran is 1-based; we’ll store iy as 1-based by padding index 0
	iy = [0] * (len(genome) + 1)
	for i, ch in enumerate(genome, start=1):
		try:
			iy[i] = _BASE_TO_IY[ch]
		except KeyError:
			raise ValueError(f"problem with genome file at position {i}: {ch!r}")
	return iy


def codon_index_from_iy(k1: int, k2: int, k3: int) -> int:
	"""
	Matches Fortran:
	  lx = sum_{ir=1..3} mm**(llr-ir)*k(ir)
	with mm=4, llr=3, k(ir) in {1..4}
	=> lx = 4**2*k1 + 4**1*k2 + 4**0*k3
	This produces indices in [21..84] (as used in the code), not [0..63].
	"""
	return (16 * k1) + (4 * k2) + k3


# -------------------------
# Core JSCB computations
# -------------------------

def compute_gene_codonfreq_and_aminofreq(
	genes: List[Gene],
	iy: List[int],
	iaa: List[List[int]],
	iab: List[int],
	ih0: int = 21,
	ih: int = 84,
) -> Tuple[List[List[int]], List[List[int]], List[int]]:
	"""
	Returns:
	  ncodonfreq[geneIndex][codonIndex] (1-based gene index; codonIndex 0..ih)
	  naminofreq[geneIndex][codonIndex]
	  gene_len_internal[geneIndex] = irend - lend - 2 (Fortran's len(j))
	"""
	ngenea = len(genes)
	# allocate with padding for 1-based gene indexing and codon index up to ih
	ncodonfreq = [[0] * (ih + 1) for _ in range(ngenea + 1)]
	naminofreq = [[0] * (ih + 1) for _ in range(ngenea + 1)]
	gene_len_internal = [0] * (ngenea + 1)

	# init
	for g in range(1, ngenea + 1):
		for c in range(ih0, ih + 1):
			ncodonfreq[g][c] = 1   # pseudocount like Fortran
			naminofreq[g][c] = 0

	for g, gene in enumerate(genes, start=1):
		lend = gene.start
		irend = gene.end
		gene_len_internal[g] = irend - lend - 2

		ncodons = (irend - lend - 2) // 3
		if gene.strand == '+':
			for ii in range(1, ncodons + 1):
				# Fortran: k(ir)=iy(lend + ii*3 - (3-ir) - 1)
				# => positions: lend+ii*3-3, -2, -1
				k1 = iy[lend + ii * 3 - 3]
				k2 = iy[lend + ii * 3 - 2]
				k3 = iy[lend + ii * 3 - 1]
				lx = codon_index_from_iy(k1, k2, k3)
				ncodonfreq[g][lx] += 1

		elif gene.strand == '-':
			for ii in range(1, ncodons + 1):
				# Fortran reads bases from the end, complements, then forms lx
				# iy(irend - ii*3 + (3-ir) + 1) for ir=1..3:
				# ir=1 => irend - 3*ii + 3
				# ir=2 => irend - 3*ii + 2
				# ir=3 => irend - 3*ii + 1
				p1 = iy[irend - 3 * ii + 3]
				p2 = iy[irend - 3 * ii + 2]
				p3 = iy[irend - 3 * ii + 1]
				k1 = _COMP_IY[p1]
				k2 = _COMP_IY[p2]
				k3 = _COMP_IY[p3]
				lx = codon_index_from_iy(k1, k2, k3)
				ncodonfreq[g][lx] += 1
		else:
			raise ValueError(f"Unexpected strand: {gene.strand!r}")

		# Build naminofreq by summing over synonymous codons
		for cod in range(ih0, ih + 1):
			found = False
			for aa_idx in range(20):
				for syn_i in range(iab[aa_idx]):
					if iaa[aa_idx][syn_i] == cod:
						total = 0
						for syn_j in range(iab[aa_idx]):
							cod2 = iaa[aa_idx][syn_j]
							total += ncodonfreq[g][cod2]
						naminofreq[g][cod] += total
						found = True
						break
				if found:
					break

	return ncodonfreq, naminofreq, gene_len_internal


def rentro_cluster(
	cluster_codon: List[int],
	cluster_am: List[int],
	icluslen: int,
	iaa: List[List[int]],
	iab: List[int],
) -> float:
	# Fortran rentro(id)
	val = 0.0
	for aa_idx in range(20):
		rentroi = 0.0
		for syn_i in range(iab[aa_idx]):
			cod = iaa[aa_idx][syn_i]
			pp = cluster_codon[cod] / float(cluster_am[cod])
			rentroi += pp * math.log(pp)
		rentroi *= (cluster_am[iaa[aa_idx][0]] / float(icluslen))
		val += rentroi
	return -val


def rentro_gene(
	gene_codon: List[int],
	gene_am: List[int],
	gene_len_codons: int,
	iaa: List[List[int]],
	iab: List[int],
) -> float:
	# Fortran rentro2(id)
	val = 0.0
	for aa_idx in range(20):
		rentroi = 0.0
		for syn_i in range(iab[aa_idx]):
			cod = iaa[aa_idx][syn_i]
			pp = gene_codon[cod] / float(gene_am[cod])
			rentroi += pp * math.log(pp)
		rentroi *= (gene_am[iaa[aa_idx][0]] / float(gene_len_codons))
		val += rentroi
	return -val


def rentrotot_cluster_cluster(
	c1_codon: List[int],
	c1_am: List[int],
	c2_codon: List[int],
	c2_am: List[int],
	lengtot: int,
	iaa: List[List[int]],
	iab: List[int],
) -> float:
	# Fortran rentrotot(id1,id2,lengtot)
	val = 0.0
	for aa_idx in range(20):
		rentroi = 0.0
		denom = float(c1_am[iaa[aa_idx][0]] + c2_am[iaa[aa_idx][0]])
		for syn_i in range(iab[aa_idx]):
			cod = iaa[aa_idx][syn_i]
			pp = (c1_codon[cod] + c2_codon[cod]) / denom
			rentroi += pp * math.log(pp)
		pp0 = (c1_am[iaa[aa_idx][0]] + c2_am[iaa[aa_idx][0]]) / float(lengtot)
		rentroi *= pp0
		val += rentroi
	return -val


def rentrotot_gene_cluster(
	gene_codon: List[int],
	gene_am: List[int],
	cl_codon: List[int],
	cl_am: List[int],
	lengtot: int,
	iaa: List[List[int]],
	iab: List[int],
) -> float:
	# Fortran rentrotot2(id1,id2,lengtot) (id1 is gene, id2 is cluster)
	val = 0.0
	for aa_idx in range(20):
		rentroi = 0.0
		denom = float(gene_am[iaa[aa_idx][0]] + cl_am[iaa[aa_idx][0]])
		for syn_i in range(iab[aa_idx]):
			cod = iaa[aa_idx][syn_i]
			pp = (gene_codon[cod] + cl_codon[cod]) / denom
			rentroi += pp * math.log(pp)
		pp0 = (gene_am[iaa[aa_idx][0]] + cl_am[iaa[aa_idx][0]]) / float(lengtot)
		rentroi *= pp0
		val += rentroi
	return -val


def divgn_cluster_cluster(
	c1_codon: List[int], c1_am: List[int], icluslen1: int,
	c2_codon: List[int], c2_am: List[int], icluslen2: int,
	lengtot: int,
	iaa: List[List[int]], iab: List[int],
) -> float:
	# Fortran divgn(i,j,lengtot)
	hhh = 2.0
	xxx = 1.0 / math.log(hhh)
	return xxx * (
		rentrotot_cluster_cluster(c1_codon, c1_am, c2_codon, c2_am, lengtot, iaa, iab)
		- (icluslen1 * rentro_cluster(c1_codon, c1_am, icluslen1, iaa, iab)
		   + icluslen2 * rentro_cluster(c2_codon, c2_am, icluslen2, iaa, iab)) / float(lengtot)
	)


def divgn2_gene_cluster(
	gene_codon: List[int], gene_am: List[int], gene_len_codons: int,
	cl_codon: List[int], cl_am: List[int], icluslen: int,
	lengtot: int,
	iaa: List[List[int]], iab: List[int],
) -> float:
	# Fortran divgn2(i,j,lengtot)
	hhh = 2.0
	xxx = 1.0 / math.log(hhh)
	return xxx * (
		rentrotot_gene_cluster(gene_codon, gene_am, cl_codon, cl_am, lengtot, iaa, iab)
		- ((gene_len_codons) * rentro_gene(gene_codon, gene_am, gene_len_codons, iaa, iab)
		   + icluslen * rentro_cluster(cl_codon, cl_am, icluslen, iaa, iab)) / float(lengtot)
	)


def cluster_jscb(
	ncodonfreq: List[List[int]],
	naminofreq: List[List[int]],
	gene_len_internal: List[int],
	iaa: List[List[int]],
	iab: List[int],
	ih0: int = 21,
	ih: int = 84,
	thres: float = 0.8,
	thres2: float = 0.0,
	thres3: float = 0.995,
	clusthres: float = 1.0,
) -> Tuple[List[int], List[List[int]], List[List[int]], List[int], List[int]]:
	"""
	Returns:
	  gene_to_cluster: 1-based gene index -> cluster id (1..nclus)
	  clusters: list of clusters; each is list of gene indices (1-based)
	  cluster_codonfreq: cluster id -> codon counts
	  cluster_aminofreq: cluster id -> amino-grouped counts
	  cluster_len_codons: cluster id -> total length in codons (icluslen)
	"""
	ngenea = len(gene_len_internal) - 1

	# --- sequential build
	clusters: List[List[int]] = []
	cluster_codonfreq: List[List[int]] = [[0] * (ih + 1)]  # dummy at index 0
	cluster_aminofreq: List[List[int]] = [[0] * (ih + 1)]
	cluster_len_codons: List[int] = [0]

	def recompute_cluster_len(cid: int) -> int:
		# Fortran: sum over 20 amino acids of nclusamfreq(cid, iaa(ii,1))
		s = 0
		am = cluster_aminofreq[cid]
		for aa_idx in range(20):
			s += am[iaa[aa_idx][0]]
		return s

	def init_cluster_with_gene(gene_idx: int) -> None:
		clusters.append([gene_idx])
		cluster_codonfreq.append([0] * (ih + 1))
		cluster_aminofreq.append([0] * (ih + 1))
		cid = len(clusters)
		for c in range(ih0, ih + 1):
			cluster_codonfreq[cid][c] = ncodonfreq[gene_idx][c]
			cluster_aminofreq[cid][c] = naminofreq[gene_idx][c]
		cluster_len_codons.append(recompute_cluster_len(cid))

	def add_gene_to_cluster(cid: int, gene_idx: int) -> None:
		clusters[cid - 1].append(gene_idx)
		for c in range(ih0, ih + 1):
			cluster_codonfreq[cid][c] += ncodonfreq[gene_idx][c]
			cluster_aminofreq[cid][c] += naminofreq[gene_idx][c]
		cluster_len_codons[cid] = recompute_cluster_len(cid)

	init_cluster_with_gene(1)
	current_cid = 1

	for gene_idx in range(2, ngenea + 1):
		gene_len_codons = gene_len_internal[gene_idx] // 3
		lengtot = cluster_len_codons[current_cid] + gene_len_codons
		ddiv = divgn2_gene_cluster(
			ncodonfreq[gene_idx], naminofreq[gene_idx], gene_len_codons,
			cluster_codonfreq[current_cid], cluster_aminofreq[current_cid], cluster_len_codons[current_cid],
			lengtot, iaa, iab
		)
		pbc = sig(float(lengtot), abs(ddiv))
		if pbc < thres:
			add_gene_to_cluster(current_cid, gene_idx)
		else:
			init_cluster_with_gene(gene_idx)
			current_cid = len(clusters)

	# --- merge by composition (thres3)
	changed = True
	while changed:
		changed = False
		nclus = len(clusters)
		outer_break = False
		for i in range(1, nclus):
			for j in range(i + 1, nclus + 1):
				lengtot = cluster_len_codons[i] + cluster_len_codons[j]
				djs = divgn_cluster_cluster(
					cluster_codonfreq[i], cluster_aminofreq[i], cluster_len_codons[i],
					cluster_codonfreq[j], cluster_aminofreq[j], cluster_len_codons[j],
					lengtot, iaa, iab
				)
				pbc = sig(float(lengtot), djs)
				if pbc < thres3:
					# merge j into i; delete j
					clusters[i - 1].extend(clusters[j - 1])
					for c in range(ih0, ih + 1):
						cluster_codonfreq[i][c] += cluster_codonfreq[j][c]
						cluster_aminofreq[i][c] += cluster_aminofreq[j][c]
					cluster_len_codons[i] = recompute_cluster_len(i)

					# remove j structures (1-based lists with dummy)
					del clusters[j - 1]
					del cluster_codonfreq[j]
					del cluster_aminofreq[j]
					del cluster_len_codons[j]

					changed = True
					outer_break = True
					break
			if outer_break:
				break

	# --- positional refinement (operon adjacency)
	def compute_rpos() -> List[List[float]]:
		nclus = len(clusters)
		# map gene -> cluster
		gene_to_cluster_local = [0] * (ngenea + 1)
		for cid, genes in enumerate(clusters, start=1):
			for g in genes:
				gene_to_cluster_local[g] = cid

		rpos = [[0.0] * (nclus + 1) for _ in range(nclus + 1)]
		for cid, genes in enumerate(clusters, start=1):
			for g in genes:
				left = g - 1
				right = g + 1
				if left < 1 or right > ngenea:
					continue
				c_left = gene_to_cluster_local[left]
				c_right = gene_to_cluster_local[right]
				if c_left == c_right and c_left != 0:
					rpos[cid][c_left] += 1.0
		return rpos

	while True:
		nclus = len(clusters)
		rpos = compute_rpos()
		merged_any = False
		for i in range(1, nclus):
			for j in range(i + 1, nclus + 1):
				rp = (rpos[i][j] / len(clusters[i - 1]) + rpos[j][i] / len(clusters[j - 1])) / 2.0
				if rp > clusthres:
					# merge j into i
					clusters[i - 1].extend(clusters[j - 1])
					for c in range(ih0, ih + 1):
						cluster_codonfreq[i][c] += cluster_codonfreq[j][c]
						cluster_aminofreq[i][c] += cluster_aminofreq[j][c]
					cluster_len_codons[i] = recompute_cluster_len(i)

					del clusters[j - 1]
					del cluster_codonfreq[j]
					del cluster_aminofreq[j]
					del cluster_len_codons[j]

					merged_any = True
					break
			if merged_any:
				break
		if not merged_any:
			break

	# --- neighborhood-based reassignment (thres2 == 0.0 means no moves in practice)
	# Build gene -> cluster lookup
	def build_gene_to_cluster() -> List[int]:
		g2c = [0] * (ngenea + 1)
		for cid, genes in enumerate(clusters, start=1):
			for g in genes:
				g2c[g] = cid
		return g2c

	gene_to_cluster = build_gene_to_cluster()

	# mimic the Fortran pass: for each cluster i, examine each gene g
	for cid in range(1, len(clusters) + 1):
		genes = list(clusters[cid - 1])
		moved_out: set[int] = set()

		for g in genes:
			if g == 1 or g == ngenea:
				continue
			c_left = gene_to_cluster[g - 1]
			c_right = gene_to_cluster[g + 1]
			if c_left == c_right and c_left != cid:
				target = c_left
				gene_len_codons = gene_len_internal[g] // 3
				lengtot = gene_len_codons + cluster_len_codons[target]
				gdjs = divgn2_gene_cluster(
					ncodonfreq[g], naminofreq[g], gene_len_codons,
					cluster_codonfreq[target], cluster_aminofreq[target], cluster_len_codons[target],
					lengtot, iaa, iab
				)
				pbc = sig(float(lengtot), gdjs)
				if pbc < thres2:
					# move gene g from cid -> target
					clusters[target - 1].append(g)
					for c in range(ih0, ih + 1):
						cluster_codonfreq[target][c] += ncodonfreq[g][c]
						cluster_aminofreq[target][c] += naminofreq[g][c]
					cluster_len_codons[target] = recompute_cluster_len(target)
					moved_out.add(g)

		if moved_out:
			clusters[cid - 1] = [g for g in clusters[cid - 1] if g not in moved_out]
			# recompute cid counts from scratch (simpler + safe)
			cluster_codonfreq[cid] = [0] * (ih + 1)
			cluster_aminofreq[cid] = [0] * (ih + 1)
			for g in clusters[cid - 1]:
				for c in range(ih0, ih + 1):
					cluster_codonfreq[cid][c] += ncodonfreq[g][c]
					cluster_aminofreq[cid][c] += naminofreq[g][c]
			cluster_len_codons[cid] = recompute_cluster_len(cid)

			gene_to_cluster = build_gene_to_cluster()

	gene_to_cluster = build_gene_to_cluster()
	return gene_to_cluster, clusters, cluster_codonfreq, cluster_aminofreq, cluster_len_codons


# -------------------------
# End-to-end pipeline (GenBank -> clusters -> GIs)
# -------------------------

def sanitize_genome(seq: str) -> str:
	nucs = ['A', 'T', 'C', 'G']
	out = []
	for ch in seq.upper():
		if ch in _BASE_TO_IY:
			out.append(ch)
		else:
			out.append(random.choice(nucs))
	return "".join(out)


def load_from_genbank(path: str) -> Tuple[str, List[Gene], List[Tuple[int, int]]]:
	genes: List[Gene] = []
	rrna: List[Tuple[int, int]] = []

	recs = list(SeqIO.parse(path, "genbank"))
	if not recs:
		raise ValueError("No GenBank records found")

	# This repo implicitly assumes one record; match that behavior.
	rec = recs[0]
	genome = sanitize_genome(str(rec.seq))

	for feature in rec.features:
		start = int(feature.location.start) + 1
		end = int(feature.location.end)
		s = feature.location.strand
		if s == 1:
			strand = "+"
		elif s == -1:
			strand = "-"
		else:
			# If Prokka omits strand (rare for CDS), default to '+'
			strand = "+"
		if feature.type == "CDS":
			if "join" not in str(feature.location):
				genes.append(Gene(strand=strand, start=start, end=end))
		elif feature.type == "rRNA":
			rrna.append((start, end))

	# match jscb.py: keep only CDS with length % 3 == 0
	filtered_genes = [g for g in genes if (g.length_nt % 3) == 0]
	return genome, filtered_genes, rrna


def write_output_gi(
	out_path: str,
	gene_to_cluster: List[int],
	cluster_sizes: Dict[int, int],
	gene_lena: List[int],
) -> None:
	# Roughly matches the info written by Fortran to JSCB_output.gi:
	# write(41,*) m, i, icsize(i), lena(m)
	# Here m is original gene index. We output:
	# geneIndex clusterId clusterSize geneLength
	with open(out_path, "w", encoding="utf-8") as f:
		f.write(f"gene_number\tclus_id\tclus_size\tgene_length\n")
		for g in range(1, len(gene_to_cluster)):
			cid = gene_to_cluster[g]
			csize = cluster_sizes.get(cid, 0) if cid > 0 else 0
			f.write(f"{g}\t{cid}\t{csize}\t{gene_lena[g]}\n")


def call_genomic_islands_like_jscb_py(
	out_tsv: str,
	gene_to_cluster: List[int],
	genes: List[Gene],
	rrna: List[Tuple[int, int]],
	min_genes: int = 8,
) -> None:
	# Find native cluster = max cluster size (by gene count)
	counts: Dict[int, int] = {}
	for g in range(1, len(gene_to_cluster)):
		cid = gene_to_cluster[g]
		if cid > 0:
			counts[cid] = counts.get(cid, 0) + 1
	if not counts:
		with open(out_tsv, "w", encoding="utf-8") as outfile:
			outfile.write("No Genomic Islands Identified\n")
		return
	native_cluster = max(counts, key=counts.get)

	actstart = [0] + [g.start for g in genes]  # 1-based
	actend = [0] + [g.end for g in genes]	  # 1-based

	prv = -1
	gmode = 0
	gicounter = 0

	with open(out_tsv, "w", encoding="utf-8") as outfile:
		for gene_idx in range(1, len(gene_to_cluster)):
			cid = gene_to_cluster[gene_idx]
			if cid != native_cluster:
				if gene_idx != prv + 1:
					startgeneno = gene_idx
					startcoord = actstart[gene_idx]
					gmode = 1
				prv = gene_idx
			elif gmode == 1:
				endgeneno = gene_idx
				gilength = endgeneno - startgeneno
				endcoord = actend[gene_idx - 1]

				if gilength >= min_genes:
					rnais = 0
					for x, y in rrna:
						# exact same containment test as jscb.py
						if startcoord < x and startcoord < y and endcoord > x and endcoord > y:
							rnais = 1
							break
					if rnais == 0:
						gicounter += 1
						if gicounter == 1:
							outfile.write("clus_id\tstart\tend\n")
						outfile.write(f"GI-{gicounter}\t{startcoord}\t{actend[gene_idx - 1]}\n")

				gmode = 0

		if gicounter == 0:
			outfile.write("No Genomic Islands Identified\n")


def main(argv: List[str]) -> int:
	args = argv[1:]
	debug = False
	if "--debug" in args:
		debug = True
		args = [a for a in args if a != "--debug"]

	if len(args) < 1:
		print("Error: Please provide input file in genbank format!!", file=sys.stderr)
		return 2

	gb_path = args[0]
	print(gb_path)
	script_path = os.path.abspath(os.path.dirname(__file__))
	rRNA_iterator = SeqIO.parse(os.path.join(script_path, 'E_coli_16S.gb'), 'genbank')
	rRNA_record = next(rRNA_iterator)
	combined_record = None
	CDS_locus_tags = list() # Ordered list of CDS locus_tags so that we can assign them back to gene numbers later
	for seq_record in SeqIO.parse(gb_path, 'genbank'):
		for feature in seq_record.features:
			if feature.type == 'CDS':
				locus_tag = feature.qualifiers['locus_tag'][0]
				CDS_locus_tags.append(locus_tag)
		if combined_record is None:
			combined_record = seq_record
		else:
			combined_record = combined_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + rRNA_record + seq_record
	combined_gbk_path = 'combined.gbk'
	SeqIO.write(combined_record, combined_gbk_path, 'genbank')
	genome, genes_raw, rrna = load_from_genbank(combined_gbk_path)
	genes_all = genes_raw

	# The Fortran filters out short genes by:
	# lena = end - start - 2; if (lena-3 > mingenelen) keep, with mingenelen=3
	# => keep if lena > 6, i.e. (end-start-2) > 6 => length_nt > 8
	# We'll implement the same.
	mingenelen = 3
	genes_used: List[Gene] = []
	for g in genes_all:
		lena = g.end - g.start - 2
		if (lena - 3) > mingenelen:
			genes_used.append(g)

	if not genes_used:
		print("No genes passed filters (CDS length % 3 == 0 and minimum length).", file=sys.stderr)
		return 1

	iy = genome_to_iy(genome)
	iaa, iab = build_iaa_iab()

	ncodonfreq, naminofreq, gene_len_internal = compute_gene_codonfreq_and_aminofreq(
		genes=genes_used, iy=iy, iaa=iaa, iab=iab
	)

	gene_to_cluster_used, clusters, _, _, _ = cluster_jscb(
		ncodonfreq=ncodonfreq,
		naminofreq=naminofreq,
		gene_len_internal=gene_len_internal,
		iaa=iaa,
		iab=iab,
		thres=0.8,
		thres2=0.0,
		thres3=0.995,
		clusthres=1.0,
	)

	# map clustered genes_used back to genes_all by exact coordinates+strand
	cluster_by_coord: Dict[Tuple[int, int, str], int] = {}
	for idx_used, g in enumerate(genes_used, start=1):
		cluster_by_coord[(g.start, g.end, g.strand)] = gene_to_cluster_used[idx_used]

	gene_to_cluster_all = [0] * (len(genes_all) + 1)
	for idx_all, g in enumerate(genes_all, start=1):
		gene_to_cluster_all[idx_all] = cluster_by_coord.get((g.start, g.end, g.strand), 0)

	cluster_sizes: Dict[int, int] = {}
	for idx_used in range(1, len(gene_to_cluster_used)):
		cid = gene_to_cluster_used[idx_used]
		if cid > 0:
			cluster_sizes[cid] = cluster_sizes.get(cid, 0) + 1

	# match Fortran lena(m) = end - start - 2 in genes_all index space
	gene_lena = [0] + [(g.end - g.start - 2) for g in genes_all]

	if debug:
		mapped_nonzero = sum(1 for cid in gene_to_cluster_all[1:] if cid > 0)
		unmapped = len(genes_all) - mapped_nonzero
		native_cid = 0
		native_size = 0
		for cid, csize in cluster_sizes.items():
			if csize > native_size:
				native_cid = cid
				native_size = csize
		print(
			f"DEBUG stats: genes_all={len(genes_all)} genes_used={len(genes_used)} "
			f"mapped={mapped_nonzero} unmapped={unmapped} nclusters={len(clusters)} "
			f"native_cluster={native_cid} native_size={native_size}",
			file=sys.stderr,
		)

	write_output_gi("JSCB_output.gi", gene_to_cluster_all, cluster_sizes, gene_lena)
	call_genomic_islands_like_jscb_py("JSCB_output.tsv", gene_to_cluster_all, genes_all, rrna)
	
	# Post Processing
	genomic_islands = {}
	jscb_output_path = 'JSCB_output.tsv'
	with open(jscb_output_path) as jscb_output:
		for i,line in enumerate(jscb_output):
			if i > 0:
				line_list = line.rstrip().split('\t')
				gi_name = line_list[0]
				start_coordinate = int(line_list[1])
				end_coordinate = int(line_list[2])
				gi_set = set(range(start_coordinate, end_coordinate + 1))
				genomic_islands[gi_name] = gi_set
				
	genomic_island_genes = dict() # keyed by GI name, holds sets of locus_tags
	for gi_name in genomic_islands:
		genomic_island_genes[gi_name] = set()

	for seq_record in SeqIO.parse(combined_gbk_path, 'genbank'):
		for feature in seq_record.features:
			if feature.type == 'gene':
				# Work out if the gene overlaps with any genomic islands
				gene_coordinate_set = set(range(int(feature.location.start), int(feature.location.end) + 1))
				locus_tag = feature.qualifiers['locus_tag'][0] # Here we take just the first value
				for gi_name in genomic_islands:
					if len(gene_coordinate_set.intersection(genomic_islands[gi_name])) > 0:
						genomic_island_genes[gi_name].add(locus_tag)

	# Now we go through the original genbank file, noting the coordinates of genes and which contigs they
	# belong to. Note, we do not assume that genomic islands won't occur at contig boundaries because
	# the JSCB program did not see where these boundaries are.
	genomic_island_coordinates = dict() # Keyed by gi_name, then contig, holds two member list [lowest_coord, highest_coord] derived from component genes
	genomic_island_genes_by_contig = dict() # Keyed by gi_name, then contig, holds lists of locus_tags

	for gi_name in genomic_islands:
		genomic_island_coordinates[gi_name] = dict()
		genomic_island_genes_by_contig[gi_name] = dict()

	for seq_record in SeqIO.parse(gb_path, 'genbank'):
		contig_name = seq_record.id
		for feature in seq_record.features:
			if feature.type == 'gene':
				locus_tag = feature.qualifiers['locus_tag'][0]
				for gi_name in genomic_island_genes:
					if locus_tag in genomic_island_genes[gi_name]:
						start_coordinate = int(feature.location.start)
						end_coordinate = int(feature.location.end)
						if contig_name in genomic_island_coordinates[gi_name]:
							if start_coordinate < genomic_island_coordinates[gi_name][contig_name][0]:
								genomic_island_coordinates[gi_name][contig_name][0] = start_coordinate
							if end_coordinate > genomic_island_coordinates[gi_name][contig_name][1]:
								genomic_island_coordinates[gi_name][contig_name][1] = end_coordinate
						else:
							genomic_island_coordinates[gi_name][contig_name] = [ start_coordinate, end_coordinate ]

						if contig_name in genomic_island_genes_by_contig[gi_name]:
							genomic_island_genes_by_contig[gi_name][contig_name].append(locus_tag)
						else:
							genomic_island_genes_by_contig[gi_name][contig_name] = [ locus_tag ]

	# Now we make a human-readable output table
	output_table_path ='genomic_islands_summary.tsv'
	output_table = open(output_table_path, 'w')
	output_table.write('GI_ID\tcontig\tstart\tend\tgenes\n')

	for gi_name in genomic_island_coordinates:
		if len(genomic_island_coordinates[gi_name]) > 1:
			# This means the GI spans multiple contigs
			# Because the order of contigs in the original file is likely random,
			# This means in reality there are two GIs
			GI_counter = 0
			for contig_name in genomic_island_coordinates[gi_name]:
				GI_counter += 1
				output_gi_name = gi_name + '_' + str(GI_counter)
				gene_string = ','.join(genomic_island_genes_by_contig[gi_name][contig_name])
				output_string = '\t'.join([ output_gi_name, contig_name, str(genomic_island_coordinates[gi_name][contig_name][0]), str(genomic_island_coordinates[gi_name][contig_name][1]), gene_string ])
				output_table.write(output_string + '\n')
		else:
			for contig_name in genomic_island_coordinates[gi_name]:
				gene_string = ','.join(genomic_island_genes_by_contig[gi_name][contig_name])
				output_string = '\t'.join([ gi_name, contig_name, str(genomic_island_coordinates[gi_name][contig_name][0]), str(genomic_island_coordinates[gi_name][contig_name][1]), gene_string ])
				output_table.write(output_string + '\n')
	output_table.close()

	# We also make a human readable version of JSCB_output.gi
	output_table_path = 'hgt_genes_summary.tsv'
	output_table = open(output_table_path, 'w')
	output_table.write('locus_tag\tcluster_id\tcluster_size\tgene_length\n')
	jscb_output_path = 'JSCB_output.gi'

	with open(jscb_output_path) as jscb_output:
		jscb_output.readline() # throw away header
		for line in jscb_output:
			line_list = line.rstrip().split()
			gene_index = int(line_list[0]) - 1

			locus_tag = CDS_locus_tags[gene_index]

			output_string = '\t'.join([locus_tag] + line_list[1:])
			output_table.write(output_string + '\n')

	output_table.close()


	return 0


if __name__ == "__main__":
	raise SystemExit(main(sys.argv))