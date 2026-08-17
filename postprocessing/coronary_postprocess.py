#!/usr/bin/env python3
#
# HCT-CCO post-processing.
# Copyright (C) 2024-2026 S. L. Vajire, J. S. Choy, G. S. Kassab and L.-C. Lee
# (Michigan State University; California Medical Innovations Institute).
#
# This file is part of HCT-CCO, a coronary-specific extension of OpenCCO
# (Copyright (C) 2023 B. Kerautret, Phuc Ngo, N. Passat, H. Talbot, C. Jaquet).
# It is free software: you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version. Distributed WITHOUT
# ANY WARRANTY; see the GNU General Public License <https://www.gnu.org/licenses/>.
#
"""
HCT-CCO post-processing (consolidated, single entry point).

Reads a generated tree XML (GXL from HctCco) and writes the FINAL
outputs directly, with no intermediate files:

  <prefix>_tree.vtp              tree with diameter-defined Strahler `final_order`
                                 and order-wise affine-transformed `radius`
  <prefix>_tree.obj              tube mesh of the transformed tree
  <prefix>_diameter_vs_order.png diameter vs DDS-order plot
  <prefix>_connectivity_matrix.csv / _diameter_stats.csv / _diameter_ranges.csv

Usage:  python3 coronary_postprocess.py  tree.xml  [output_prefix]

The DDS-ordering, order-wise affine diameter transform and tube-meshing are the
exact routines used in the study (merged here from the ordering, line-graph and
vtp->obj scripts). Only the front-end changed: the tree XML is parsed directly,
which removes the Swift graph-export + intermediate .dat / .vtp hand-offs.
"""
import sys, os, re
import numpy as np
import vtk
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Headless one-shot run: save figures to file instead of opening a window.
_PLOT_PATH = ["diameter_vs_order.png"]
def _save_show(*args, **kwargs):
    plt.savefig(_PLOT_PATH[0], dpi=150, bbox_inches="tight")
    plt.close("all")
plt.show = _save_show


def build_polydata_from_tree_xml(xml_path):
    """Parse a HctCco GXL tree and build a vtkPolyData of poly-lines with
    a CellData 'radius' array (one value per segment). Replaces the Swift
    graph-export + line-graph steps by reading the tree XML directly."""
    text = open(xml_path).read()

    node_pos = {}
    for m in re.finditer(r'<node id="n(\d+)">(.*?)</node>', text, re.S):
        nid = int(m.group(1))
        pm = re.search(r'name="\s*position">\s*<tup>\s*'
                       r'<float>([-\d.eE+]+)</float>\s*'
                       r'<float>([-\d.eE+]+)</float>\s*'
                       r'<float>([-\d.eE+]+)</float>', m.group(2))
        if pm:
            node_pos[nid] = (float(pm.group(1)), float(pm.group(2)), float(pm.group(3)))

    edges = []
    for m in re.finditer(r'<edge id="e\d+" to="n(\d+)" from="n(\d+)">(.*?)</edge>',
                         text, re.S):
        to, frm, body = int(m.group(1)), int(m.group(2)), m.group(3)
        rm = re.search(r'name="\s*radius">\s*<float>([-\d.eE+]+)</float>', body)
        radius = float(rm.group(1)) if rm else 0.0
        edges.append((frm, to, radius))

    max_id = max(node_pos) if node_pos else -1
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(max_id + 1)
    for nid in range(max_id + 1):
        x, y, z = node_pos.get(nid, (0.0, 0.0, 0.0))
        points.SetPoint(nid, x, y, z)

    lines = vtk.vtkCellArray()
    radius_arr = vtk.vtkDoubleArray()
    radius_arr.SetName("radius")
    radius_arr.SetNumberOfComponents(1)
    for frm, to, r in edges:
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, frm)
        line.GetPointIds().SetId(1, to)
        lines.InsertNextCell(line)
        radius_arr.InsertNextValue(r)

    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    poly.SetLines(lines)
    poly.GetCellData().AddArray(radius_arr)
    return poly


# ============================================================================
#  Study routines below are reproduced verbatim (DDS ordering, affine transform,
#  morphometry CSVs, tube meshing). Only the driver at the bottom is new.
# ============================================================================

MU = np.array([0.009, 0.013, 0.0187, 0.0346, 0.0716,
               0.150, 0.303, 0.467, 0.715, 1.492, 3.176])

SD = np.array([0.0009, 0.0017, 0.0026, 0.0074, 0.0172,
               0.035, 0.054, 0.056, 0.130, 0.300, 0.600])

def _bounds(ref_order: int):
    """Return (lo, hi) of μ±σ band for *reference* order (0-based)."""
    lo = MU[ref_order] - SD[ref_order]
    hi = MU[ref_order] + SD[ref_order]
    return lo, hi

def _build_ladders(tree_orders: np.ndarray, offset: int):
    """Pre-compute {order → [next_diam, step]} ladders inside shifted bands."""
    counts = np.bincount(tree_orders)               # orders start at 1
    ladders = {}
    for o in np.nonzero(counts)[0]:                 # o is 0-based tree order
        n = counts[o]
        lo, hi = _bounds(o + offset)
        step = 0.0 if n == 1 else (hi - lo) / (n - 1)
        ladders[o] = [hi, step]                     # pointer starts at top
    return ladders

def ladder_assign(root, children, orders, offset, eps=1e-4):
    """
    Deterministic μ±σ-bounded, monotone radii.
      root       : int                     id of root segment
      children   : list[list[int]]         children[i] = list of child ids
      orders     : np.ndarray[int]         diameter-defined order (1-based)
      offset     : int                     0-2 so root maps → ref-order 11
    Returns → np.ndarray of new **diameters** (same length as orders).
    """
    ladders = _build_ladders(orders - 1, offset)    # tree orders 0-based
    diam = np.empty_like(orders, dtype=float)       # will hold diameters

    stack = [(root, np.inf)]                        # (seg_id, parent_diam)
    while stack:
        seg, parent_d = stack.pop()
        o = orders[seg] - 1                         # tree order 0-based
        nxt, step = ladders[o]

        # highest rung that satisfies taper
        d = nxt
        while d >= parent_d - eps and step > 0.0:
            d -= step
        lo, _ = _bounds(o + offset)
        d = max(d, lo)                              # clamp safety
        ladders[o][0] = d - step                    # move pointer

        diam[seg] = d                              # store diameter
        for c in children[seg]:
            stack.append((c, d))
    return diam



def affine_rescale_by_order(diameters: np.ndarray, orders: np.ndarray, offset: int,
                            clip_min: float = 1e-12, ddof: int = 1):
    """Order-wise affine calibration to match reference mean±SD.

    For each *tree* order o (1-based):
      - map to reference order r = o + offset (also 1-based)
      - compute simulated (mu_s, sd_s) over all segments with that order
      - compute reference (mu_r, sd_r) from MU/SD
      - define b_o = sd_r / sd_s, a_o = mu_r - b_o * mu_s
      - rescale every diameter in that order: d' = a_o + b_o*d

    Returns:
      new_diameters : np.ndarray[float] same shape as diameters
      params        : dict[int, tuple[float,float]] mapping tree_order -> (a_o, b_o)
    """
    diameters = np.asarray(diameters, dtype=float)
    orders    = np.asarray(orders, dtype=int)
    new_d = diameters.copy()

    params = {}
    max_ref = len(MU)

    for o in np.unique(orders):
        if o <= 0:
            continue
        idx = np.where(orders == o)[0]
        if idx.size == 0:
            continue

        r = int(o + offset)  # 1-based reference order
        if r < 1:
            r = 1
        if r > max_ref:
            r = max_ref

        mu_r = float(MU[r-1])
        sd_r = float(SD[r-1])

        vals = diameters[idx]
        mu_s = float(np.mean(vals))
        if idx.size >= 2:
            sd_s = float(np.std(vals, ddof=ddof))
        else:
            sd_s = 0.0

        # Fallbacks for degenerate cases
        if sd_s <= 0.0 or not np.isfinite(sd_s):
            # Mean-only scaling (preserves ordering, matches mean)
            if mu_s > 0 and np.isfinite(mu_s):
                b = mu_r / mu_s
                a = 0.0
            else:
                b = 0.0
                a = mu_r
        else:
            b = sd_r / sd_s
            a = mu_r - b * mu_s

        new_d[idx] = a + b * vals
        if clip_min is not None:
            new_d[idx] = np.maximum(new_d[idx], clip_min)

        params[int(o)] = (float(a), float(b))

    return new_d, params


def enforce_tapering(root: int, children: list, diameters: np.ndarray,
                     eps: float = 1e-4, clip_min: float = 1e-12):
    """Enforce parent->child tapering: child_d < parent_d - eps along the tree.

    Note: This step may slightly change order-wise mean/SD after calibration.
    """
    d = np.asarray(diameters, dtype=float).copy()
    stack = [(root, np.inf)]
    while stack:
        seg, parent_d = stack.pop()
        if np.isfinite(parent_d):
            d[seg] = min(d[seg], parent_d - eps)
        if clip_min is not None:
            d[seg] = max(d[seg], clip_min)
        for c in children[seg]:
            stack.append((c, d[seg]))
    return d


def enforce_tapering_balanced(root: int, children: list, diameters: np.ndarray,
                              orders: np.ndarray, eps: float = 1e-4, 
                              clip_min: float = 1e-12, max_iters: int = 50):
    """
    Enforce tapering while trying to preserve order-wise statistics.
    
    Strategy: When child > parent, split the correction between them
    (expand parent slightly, shrink child slightly) to minimize impact.
    """
    d = np.asarray(diameters, dtype=float).copy()
    
    # Build parent lookup
    parent_of = [-1] * len(d)
    stack = [root]
    while stack:
        seg = stack.pop()
        for c in children[seg]:
            parent_of[c] = seg
            stack.append(c)
    
    for iteration in range(max_iters):
        violations = 0
        
        # Process in BFS order (root to leaves)
        visited = [False] * len(d)
        queue = [root]
        
        while queue:
            seg = queue.pop(0)
            if visited[seg]:
                continue
            visited[seg] = True
            
            p = parent_of[seg]
            if p >= 0 and d[seg] >= d[p] - eps:
                # Violation! Split the correction
                gap_needed = eps
                current_gap = d[p] - d[seg]
                correction = (gap_needed - current_gap) / 2.0 + 1e-6
                
                # Expand parent, shrink child (balanced)
                d[p] += correction
                d[seg] -= correction
                violations += 1
            
            for c in children[seg]:
                queue.append(c)
        
        if violations == 0:
            break
    
    # Final clamp
    if clip_min is not None:
        d = np.maximum(d, clip_min)
    
    return d

def load_vascular_network(file_path):
    """
    Loads a vascular network from a VTK file.

    Parameters:
    - file_path: Path to the VTK file.

    Returns:
    - vtkPolyData object representing the vascular network.
    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file_path)
    reader.Update()
    return reader.GetOutput()

def get_vessel_diameters(vtk_polydata):
    """
    Extracts vessel diameters from the VTK PolyData.

    Parameters:
    - vtk_polydata: The vascular network data.

    Returns:
    - Numpy array of diameters for each segment.
    """
    radii = vtk.vtkDoubleArray.SafeDownCast(vtk_polydata.GetCellData().GetArray("radius"))
    if radii is None:
        raise ValueError("The VTK file does not contain a 'radius' array in CellData.")
    return np.array([2 * radii.GetValue(i) for i in range(vtk_polydata.GetNumberOfCells())])

def build_segment_connections(vtk_polydata):
    """
    Builds parent and child relationships between segments using ONLY endpoints
    (prevents fake junctions from shared interior points).
    """
    num_cells = vtk_polydata.GetNumberOfCells()
    segment_connections = {
        'parent': [[] for _ in range(num_cells)],
        'children': [[] for _ in range(num_cells)]
    }

    # Build mapping from endpoint point IDs to segment IDs
    point_to_segment = {}
    for seg_id in range(num_cells):
        cell = vtk_polydata.GetCell(seg_id)
        if cell.GetNumberOfPoints() < 2:
            continue

        start_pid = cell.GetPointId(0)
        end_pid = cell.GetPointId(cell.GetNumberOfPoints() - 1)

        for pid in (start_pid, end_pid):
            if pid not in point_to_segment:
                point_to_segment[pid] = []
            point_to_segment[pid].append(seg_id)

    # Determine parent and children segments using the endpoint map
    for seg_id in range(num_cells):
        cell = vtk_polydata.GetCell(seg_id)
        if cell.GetNumberOfPoints() < 2:
            continue

        start_pid = cell.GetPointId(0)
        end_pid = cell.GetPointId(cell.GetNumberOfPoints() - 1)

        downstream_segments = [s for s in point_to_segment.get(end_pid, []) if s != seg_id]
        upstream_segments   = [s for s in point_to_segment.get(start_pid, []) if s != seg_id]

        segment_connections['children'][seg_id].extend(downstream_segments)
        segment_connections['parent'][seg_id].extend(upstream_segments)

    return segment_connections

def get_segment_endpoints(vtk_polydata, seg_id):
    cell = vtk_polydata.GetCell(seg_id)
    return cell.GetPointId(0), cell.GetPointId(cell.GetNumberOfPoints() - 1)

def build_point_degree_map(vtk_polydata):
    degree = {}
    for seg_id in range(vtk_polydata.GetNumberOfCells()):
        a, b = get_segment_endpoints(vtk_polydata, seg_id)
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1
    return degree                        # point-id ➜ #incident segments

def neighbor_segments(seg_id, segment_connections):
    """Return every upstream + downstream neighbour of seg_id."""
    return (
        segment_connections["parent"][seg_id]
        + segment_connections["children"][seg_id]
    )

def collapse_straight_chains(vtk_polydata, segment_connections):
    """
    Merge every maximal run of degree-2 points into one element.
    Returns
    -------
    seg2elem    : list[int]           # len = #segments
    elem_chains : list[list[int]]     # segment-ID chains
    elem_diam   : np.ndarray[float]   # representative diameter per element
    """
    # --- helper ------------------------------------------------------------
    def get_endpoints(seg_id):
        cell = vtk_polydata.GetCell(seg_id)
        return cell.GetPointId(0), cell.GetPointId(cell.GetNumberOfPoints() - 1)

    # --- point-to-segments map & degree ------------------------------------
    n_seg = vtk_polydata.GetNumberOfCells()
    p2s   = {}            # point-ID ➜ list of seg-IDs
    degree = {}
    for s in range(n_seg):
        a, b = get_endpoints(s)
        for p in (a, b):
            p2s.setdefault(p, []).append(s)
            degree[p] = degree.get(p, 0) + 1     # add one incident segment

    # --- visit flags --------------------------------------------------------
    seg2elem   = [-1] * n_seg
    elem_chains, elem_diam = [], []
    elem_id = 0

    radius_arr = vtk.vtkDoubleArray.SafeDownCast(
        vtk_polydata.GetCellData().GetArray("radius")
    )

    # --- start a walk from every branch OR leaf point ----------------------
    branch_or_leaf = [p for p, d in degree.items() if d != 2]
    for bp in branch_or_leaf:
        for seg in p2s[bp]:
            if seg2elem[seg] != -1:
                continue                        # already part of an element
            chain = [seg]
            seg2elem[seg] = elem_id
            # determine forward direction
            cur_pt  = bp
            cur_seg = seg
            nxt_pt  = get_endpoints(cur_seg)[1] if cur_pt == get_endpoints(cur_seg)[0] \
                                                else get_endpoints(cur_seg)[0]
            # walk until the next branch/leaf
            while degree[nxt_pt] == 2:
                # exactly one *other* segment shares nxt_pt
                nxt_seg = [s for s in p2s[nxt_pt] if s != cur_seg][0]
                if seg2elem[nxt_seg] != -1:
                    break                       # safety: already assigned
                chain.append(nxt_seg)
                seg2elem[nxt_seg] = elem_id
                cur_seg = nxt_seg
                a, b = get_endpoints(cur_seg)
                cur_pt, nxt_pt = nxt_pt, (b if nxt_pt == a else a)

            # store element
            elem_chains.append(chain)
            elem_diam.append(
                np.mean([2.0 * radius_arr.GetValue(s) for s in chain])
            )
            elem_id += 1

    return seg2elem, elem_chains, np.asarray(elem_diam, float)

def build_element_connections(elem_chains, seg2elem, segment_connections):
    """
    Build a *directed* acyclic element tree with genuine leaves.
    """
    Ne = len(elem_chains)
    elem_conn = {'parent': [[] for _ in range(Ne)],
                 'children': [[] for _ in range(Ne)]}

    # link every upstream segment to its downstream neighbour
    for seg_id, kids in enumerate(segment_connections['children']):
        e_parent = seg2elem[seg_id]
        for k in kids:
            e_child = seg2elem[k]
            if e_child != e_parent:
                elem_conn['children'][e_parent].append(e_child)
                elem_conn['parent'][e_child].append(e_parent)

    # drop duplicates
    elem_conn['parent']   = [list(dict.fromkeys(p)) for p in elem_conn['parent']]
    elem_conn['children'] = [list(dict.fromkeys(c)) for c in elem_conn['children']]
    return elem_conn

def propagate_to_segments(elem_orders, seg2elem, n_segments):
    SO_seg = np.zeros(n_segments, dtype=int)
    for seg_id, elem_id in enumerate(seg2elem):
        SO_seg[seg_id] = elem_orders[elem_id]
    return SO_seg

def compute_initial_strahler_orders(segment_connections, num_segments):
    """
    Computes the initial Strahler orders based on network topology.

    Parameters:
    - segment_connections: Dictionary with 'parent' and 'children' lists.
    - num_segments: Total number of segments.

    Returns:
    - Numpy array of initial Strahler orders for each segment.
    """
    SO = np.zeros(num_segments, dtype=int)
    assigned = np.zeros(num_segments, dtype=bool)
    # Start by identifying terminal segments (segments with no children)
    terminal_segments = [i for i in range(num_segments) if len(segment_connections['children'][i]) == 0]
    SO[terminal_segments] = 1
    assigned[terminal_segments] = True
    count = len(terminal_segments)

    # Iteratively assign orders to parent segments
    while count < num_segments:
        for seg_id in range(num_segments):
            if assigned[seg_id]:
                continue
            # Check if all children have assigned orders
            child_ids = segment_connections['children'][seg_id]
            if all(assigned[child_id] for child_id in child_ids):
                # All children have assigned orders
                child_orders = [SO[child_id] for child_id in child_ids]
                if len(child_orders) == 0:
                    SO[seg_id] = 1
                elif all(order == child_orders[0] for order in child_orders):
                    SO[seg_id] = child_orders[0] + 1
                else:
                    SO[seg_id] = max(child_orders)
                assigned[seg_id] = True
                count += 1

    return SO

def compute_diameter_defined_strahler_orders(SO, diameters, max_iterations=100, convergence_threshold=0.01):
    """
    Computes the diameter-defined Strahler orders.

    Parameters:
    - SO: Initial Strahler orders.
    - diameters: Diameters of segments.
    - max_iterations: Maximum number of iterations.
    - convergence_threshold: Threshold for convergence.

    Returns:
    - SOF: Final diameter-defined Strahler orders.
    - DOF: List of diameters for each order.
    - iter_count: Number of iterations performed.
    """
    Ns = len(SO)
    SOdd = SO.copy()
    MO = SO.max()
    iter_count = 0
    change = 1.0

    while change > convergence_threshold and iter_count < max_iterations:
        iter_count += 1
        Dmean = np.zeros(MO)
        Dsd = np.zeros(MO)
        DO = [ [] for _ in range(MO) ]
        # Calculate mean and std of diameters for each order
        for i in range(MO):
            indices = np.where(SOdd == (i + 1))[0]
            if len(indices) > 0:
                D_i = diameters[indices]
                DO[i] = D_i
                Dmean[i] = D_i.mean()
                Dsd[i] = D_i.std()
            else:
                Dmean[i] = 0
                Dsd[i] = 0

        # Calculate diameter boundaries
        DB = (Dmean[1:] + Dsd[1:] + Dmean[:-1] - Dsd[:-1]) / 2
        nB = Dmean[-1] + Dsd[-1]
        # Reassign orders
        SOdd_new = np.zeros(Ns, dtype=int)
        for idx in range(Ns):
            D = diameters[idx]
            if D > nB:
                new_order = MO + 1
            else:
                TF = D < DB
                nO = MO - np.sum(TF)
                new_order = nO
            SOdd_new[idx] = new_order

        change = np.sum(SOdd_new != SOdd) / Ns
        SOdd = SOdd_new.copy()
        MO = SOdd.max()
        print(f"Iteration {iter_count}, change: {change}")

    # Remove empty orders if any
    unique_orders = np.unique(SOdd)
    order_mapping = {old_order: new_order for new_order, old_order in enumerate(sorted(unique_orders), start=1)}
    SOF = np.array([order_mapping[order] for order in SOdd])
    DOF = [[] for _ in range(SOF.max())]
    for idx in range(Ns):
        DOF[SOF[idx] - 1].append(diameters[idx])
    DOF = [np.array(diameters) for diameters in DOF]

    iterF = iter_count
    return SOF, DOF, iterF

def plot_diameters(vtk_polydata, orders, title="Diameter vs Diameter-Defined Strahler Order"):
    """
    Plots diameters against their diameter-defined Strahler orders.

    Parameters:
    - vtk_polydata: The vascular network data.
    - orders: Diameter-defined Strahler orders for each segment.
    - title: Title of the plot.
    """
    diameters = []
    strahler_orders = []
    for i in range(vtk_polydata.GetNumberOfCells()):
        radius = vtk.vtkDoubleArray.SafeDownCast(vtk_polydata.GetCellData().GetArray("radius")).GetValue(i)
        diameter = 2 * radius
        diameters.append(diameter)
        strahler_orders.append(orders[i])

    plt.figure(figsize=(10, 5))
    plt.scatter(strahler_orders, diameters, alpha=0.5)
    plt.xlabel('Strahler Order')
    plt.ylabel('Diameter (units)')
    plt.title(title)
    plt.grid(True)
    plt.show()

def compute_bidirectional_connectivity_matrix(orders, vtk_polydata):
    """
    Compute *directed* parent→child connectivity between diameter-defined orders,
    using the element tree (Kassab-style connectivity matrix).

    Parameters
    ----------
    orders : np.ndarray[int]
        Diameter-defined Strahler order for each *segment* (SOF).
    vtk_polydata : vtkPolyData
        Vascular network geometry.

    Returns
    -------
    connectivity_counts : dict[(int,int) -> list[int]]
        For each pair (m, n) = (daughter_order, parent_order), a list of
        counts, one value per parent *element* of order n giving the number
        of daughter elements of order m that spring from that parent.
    """
    num_segments = vtk_polydata.GetNumberOfCells()

    # Rebuild segment-level connectivity
    segment_connections = build_segment_connections(vtk_polydata)

    # Collapse to elements exactly as in the main pipeline
    seg2elem, elem_chains, _ = collapse_straight_chains(
        vtk_polydata, segment_connections
    )

    Ne = len(elem_chains)

    # Infer element orders from segment orders: assume uniform order per element
    elem_orders = np.zeros(Ne, dtype=int)
    for e_id, chain in enumerate(elem_chains):
        if not chain:
            continue
        # all segments in a chain should share the same order
        elem_orders[e_id] = orders[chain[0]]

    # Build directed element-level connectivity
    elem_conn = build_element_connections(elem_chains, seg2elem, segment_connections)

    unique_orders = np.unique(elem_orders)
    connectivity_counts = {(m, n): [] for m in unique_orders for n in unique_orders}

    # For each parent element, count daughters of each order
    for e_parent in range(Ne):
        n = elem_orders[e_parent]
        child_ids = elem_conn["children"][e_parent]
        child_orders = [elem_orders[c] for c in child_ids]

        for m in unique_orders:
            count_m = child_orders.count(m)
            connectivity_counts[(m, n)].append(count_m)

    return connectivity_counts

def compute_mean_se(connectivity_counts):
    """
    Computes the mean and standard error for the connectivity counts.

    Parameters:
    - connectivity_counts: Dictionary with keys (m, n) and values as lists of counts.

    Returns:
    - mean_se: Dictionary with keys (m, n) and values as (mean, standard error).
    """
    mean_se = {}
    for key, counts in connectivity_counts.items():
        n_samples = len(counts)
        if n_samples > 0:
            mean = np.mean(counts)
            se = np.std(counts, ddof=1) / np.sqrt(n_samples)
        else:
            mean = 0
            se = 0
        mean_se[key] = (mean, se)
    return mean_se

def write_matrix_to_csv(mean_se, unique_orders, filename="connectivity_matrix.csv"):
    """
    Writes the (possibly non-symmetric) connectivity matrix to a CSV file.

    Parameters
    ----------
    mean_se : dict[(int,int) -> (float,float)]
        (mean, standard error) for each (m, n) = (daughter, parent).
    unique_orders : iterable[int]
        Orders present in the data.
    filename : str
        Output CSV file name.
    """
    unique_orders = sorted(unique_orders)

    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        # header
        header = ["Order m \\ Order n"] + [str(n) for n in unique_orders]
        writer.writerow(header)

        for m in unique_orders:
            row = [str(m)]
            for n in unique_orders:
                mean, se = mean_se.get((m, n), (0.0, 0.0))

                # you can keep or tweak this thresholding if you like
                if mean < 0.005:
                    cell_value = "0"
                else:
                    cell_value = f"{mean:.3f} ± {se:.3f}"

                row.append(cell_value)
            writer.writerow(row)


def compute_diameter_stats_per_order(orders, diameters):
    """
    Computes the mean and standard deviation of diameters for each order.

    Parameters:
    - orders: Diameter-defined Strahler orders for each segment.
    - diameters: Diameters of segments.

    Returns:
    - diameter_stats: Dictionary with keys as orders and values as (mean, std deviation).
    """
    unique_orders = np.unique(orders)
    diameter_stats = {}
    for order in unique_orders:
        indices = np.where(orders == order)[0]
        diameters_in_order = diameters[indices]
        mean_diameter = np.mean(diameters_in_order)
        std_diameter = np.std(diameters_in_order, ddof=1)  # Sample standard deviation
        diameter_stats[order] = (mean_diameter, std_diameter)
    return diameter_stats

def write_diameter_stats_to_csv(diameter_stats, filename="diameter_stats.csv"):
    """
    Writes the diameter statistics to a CSV file.

    Parameters:
    - diameter_stats: Dictionary with keys as orders and values as (mean, std deviation).
    - filename: Name of the CSV file.
    """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(['Order', 'Diameter (units) ± SD'])
        # Write data
        for order in sorted(diameter_stats.keys()):
            mean_diameter, std_diameter = diameter_stats[order]
            row = [order, f"{mean_diameter:.2f} ± {std_diameter:.5f}"]
            writer.writerow(row)

def compute_diameter_ranges_per_order(orders, diameters):
    """
    Computes the minimum and maximum diameters for each order.

    Parameters:
    - orders: Diameter-defined Strahler orders for each segment.
    - diameters: Diameters of segments.

    Returns:
    - diameter_ranges: Dictionary with keys as orders and values as (min_diameter, max_diameter).
    """
    unique_orders = np.unique(orders)
    diameter_ranges = {}
    for order in unique_orders:
        indices = np.where(orders == order)[0]
        diameters_in_order = diameters[indices]
        min_diameter = np.min(diameters_in_order)
        max_diameter = np.max(diameters_in_order)
        diameter_ranges[order] = (min_diameter, max_diameter)
    return diameter_ranges

def write_diameter_ranges_to_csv(diameter_ranges, filename="diameter_ranges.csv"):
    """
    Writes the diameter ranges to a CSV file.

    Parameters:
    - diameter_ranges: Dictionary with keys as orders and values as (min_diameter, max_diameter).
    - filename: Name of the CSV file.
    """
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(['Order', 'Diameter Range (units)'])
        # Write data
        for order in sorted(diameter_ranges.keys()):
            min_diameter, max_diameter = diameter_ranges[order]
            diameter_range = f"{min_diameter:.2f} - {max_diameter:.2f}"
            row = [order, diameter_range]
            writer.writerow(row)

def cap_orders_to_eleven(SOF, DOF):
    try:
        max_order = int(np.max(SOF))
    except Exception:
        return SOF, DOF
    if max_order <= 11:
        return SOF, DOF
    old_len = len(DOF)
    if old_len < max_order:
        DOF = DOF + [np.array([])] * (max_order - old_len)
    new_DOF = []
    merge_parts = []
    if old_len >= 1:
        merge_parts.append(DOF[0])
    if old_len >= 2:
        merge_parts.append(DOF[1])
    if len(merge_parts) > 0:
        has_any = any((isinstance(p, np.ndarray) and p.size > 0) for p in merge_parts)
        new_DOF.append(np.concatenate([p for p in merge_parts if isinstance(p, np.ndarray) and p.size > 0]) if has_any else np.array([]))
    else:
        new_DOF.append(np.array([]))
    for i in range(2, max_order):
        new_DOF.append(DOF[i])
    SOF_new = SOF.copy()
    SOF_new = np.asarray(SOF_new, dtype=int)
    SOF_new[SOF_new == 2] = 1
    ge3 = SOF_new >= 3
    SOF_new[ge3] = SOF_new[ge3] - 1
    return SOF_new, new_DOF

# -----------------------------------------------------------------
# Main script
# -----------------------------------------------------------------

# --- tube/obj helpers (from vtp->obj step, verbatim) ---
def ensure_point_radius(polydata: vtk.vtkPolyData, array_name: str) -> vtk.vtkPolyData:
    """If 'array_name' lives in CellData, convert it to PointData."""
    pd, cd = polydata.GetPointData(), polydata.GetCellData()
    if pd.HasArray(array_name):
        return polydata  # already in points

    if cd.HasArray(array_name):
        print("[INFO] Converting CellData → PointData for 'radius'")
        c2p = vtk.vtkCellDataToPointData()
        c2p.SetInputData(polydata)
        c2p.PassCellDataOff()
        c2p.Update()
        return c2p.GetOutput()

    print(f"[WARN] '{array_name}' missing – tubes will use fallback radius")
    return polydata


def build_tubes(polydata: vtk.vtkPolyData,
                array_name: str,
                n_sides: int,
                capping: bool,
                fallback_radius: float) -> vtk.vtkPolyData:
    """Sweep absolute-radius tubes along every poly-line."""
    use_default = not polydata.GetPointData().HasArray(array_name)

    tube = vtk.vtkTubeFilter()
    tube.SetInputData(polydata)
    tube.SetNumberOfSides(n_sides)
    tube.SetCapping(capping)

    if use_default:
        tube.SetRadius(fallback_radius)
        tube.SetVaryRadiusToOff()
    else:
        tube.SetInputArrayToProcess(
            0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, array_name
        )
        tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
        tube.SetRadius(1.0)  # ignored in absolute-scalar mode

    tube.Update()
    return tube.GetOutput()



# ============================================================================
#  Consolidated driver: tree XML -> final vtp + obj + plot + CSVs (one call).
#  Mirrors the original __main__ flow of the ordering script, but takes its
#  input from the parsed tree XML and writes the transformed .obj in the same run.
# ============================================================================
# --- hemodynamics: steady lumped-parameter (Poiseuille) check on the raw tree ---
MU_VISC_CGS  = 0.036   # blood viscosity 3.6 cP in dyn.s/cm^2
PERF_DENSITY = 3.5     # physiological (vasodilated) perfusion, mL/min/g
PERF_REST    = 1.0     # resting coronary perfusion, mL/min/g
PERF_STRESS  = 3.5     # stress (vasodilated) coronary perfusion, mL/min/g
TISSUE_DENS  = 1.05    # g/cm^3

def _read_domain_mass_g(vol_path):
    """Myocardial (mass g, volume cm3) from a DGtal .vol mask: text header ended
    by a '.' line, then one uint8 per voxel; foreground = non-zero."""
    raw = open(vol_path, "rb").read()
    he = raw.find(b"\n.\n")
    header = raw[:he].decode("latin1")
    d = {k: int(v) for k, v in re.findall(r'(X|Y|Z|Voxel-Size):\s*(-?\d+)', header)}
    X, Y, Z, vs = d['X'], d['Y'], d['Z'], d['Voxel-Size']
    data = np.frombuffer(raw[he + 3: he + 3 + X * Y * Z], dtype=np.uint8)
    V_cm3 = int((data > 0).sum()) * (vs ** 3) / 1000.0
    return V_cm3 * TISSUE_DENS, V_cm3

def _read_segments_flow_radius(xml_path):
    """(radius mm, flow) per segment from the raw tree XML, in edge order."""
    text = open(xml_path).read()
    rad, flow = [], []
    for m in re.finditer(r'<edge id="e\d+" to="n\d+" from="n\d+">(.*?)</edge>', text, re.S):
        b = m.group(1)
        r = re.search(r'name="\s*radius">\s*<float>([-\d.eE+]+)', b)
        f = re.search(r'name="\s*flow">\s*<float>([-\d.eE+]+)', b)
        rad.append(float(r.group(1)) if r else 0.0)
        flow.append(float(f.group(1)) if f else 0.0)
    return np.asarray(rad), np.asarray(flow)

def compute_hemodynamics(xml_path, prefix, domain_vol, orders):
    """Steady lumped-parameter resistive-network hemodynamics (paper Section 2.3.4).
    Each segment is a Poiseuille resistance R = 8*mu*L/(pi*r^4); each terminal is coupled to a
    distal microvascular resistance representing the unresolved downstream bed. A perfusion
    pressure of 100 mmHg is prescribed at the inlet and the total inlet flow is set to the
    physiological perfusion of the segmented mass (Q = q*rho*V); the distal resistance is
    chosen so both boundary conditions hold. The linear network is solved for per-segment flow
    and pressure, then WSS = 4*mu*Q/(pi*r^3)."""
    MMHG = 1333.22
    mass_g, V_cm3 = _read_domain_mass_g(domain_vol)
    text = open(xml_path).read()
    npos = {}
    for m in re.finditer(r'<node id="n(\d+)">(.*?)</node>', text, re.S):
        pm = re.search(r'name="\s*position">\s*<tup>\s*<float>([-\d.eE+]+)</float>\s*'
                       r'<float>([-\d.eE+]+)</float>\s*<float>([-\d.eE+]+)</float>', m.group(2))
        if pm:
            npos[int(m.group(1))] = np.array([float(pm.group(1)), float(pm.group(2)), float(pm.group(3))])
    frm, to, rad = [], [], []
    for m in re.finditer(r'<edge id="e\d+" to="n(\d+)" from="n(\d+)">(.*?)</edge>', text, re.S):
        r = re.search(r'name="\s*radius">\s*<float>([-\d.eE+]+)', m.group(3))
        to.append(int(m.group(1))); frm.append(int(m.group(2))); rad.append(float(r.group(1)))
    n = len(rad); rad = np.asarray(rad)
    if n == 0:
        print("[hemo] no segments; skipping"); return
    node2seg = {to[i]: i for i in range(n)}
    parent = [[node2seg[frm[i]]] if frm[i] in node2seg else [] for i in range(n)]
    children = [[] for _ in range(n)]
    for i in range(n):
        if parent[i]:
            children[parent[i][0]].append(i)
    root = [i for i in range(n) if not parent[i]][0]
    L = np.maximum(np.array([np.linalg.norm(npos[to[i]] - npos[frm[i]]) for i in range(n)]), 1e-3)
    Rseg = 8.0 * MU_VISC_CGS * (L / 10.0) / (np.pi * (rad / 10.0) ** 4)   # dyn.s/cm5

    Pin, Pv = 100.0, 0.0
    post = []; seen = set(); st = [(root, False)]
    while st:
        nd, dn = st.pop()
        if dn: post.append(nd); continue
        if nd in seen: continue
        seen.add(nd); st.append((nd, True))
        for c in children[nd]: st.append((c, False))

    def solve(Qtot):
        """Solve the resistive network for a prescribed total inlet flow (cm3/s)."""
        def Req_of(Rt):
            Req = np.zeros(n)
            for i in post:
                Req[i] = Rseg[i] + Rt if not children[i] else Rseg[i] + 1.0 / sum(1.0 / Req[c] for c in children[i])
            return Req
        dP = (Pin - Pv) * MMHG
        lo, hi = 0.0, 1e14                              # find distal resistance matching Pin and Qtot
        for _ in range(80):
            mid = (lo + hi) / 2.0
            if dP / Req_of(mid)[root] > Qtot: lo = mid
            else: hi = mid
        Req = Req_of((lo + hi) / 2.0)
        Q = np.zeros(n); Q[root] = Qtot                 # top-down flow split by conductance
        for i in post[::-1]:
            if children[i]:
                inv = sum(1.0 / Req[c] for c in children[i])
                for c in children[i]:
                    Q[c] = Q[i] * (1.0 / Req[c]) / inv
        P = np.zeros(n); P[root] = (Pin * MMHG - Q[root] * Rseg[root]) / MMHG   # node pressures (mmHg)
        for i in post[::-1]:
            for c in children[i]:
                P[c] = P[i] - Q[c] * Rseg[c] / MMHG
        wss = 4.0 * MU_VISC_CGS * Q / (np.pi * (rad / 10.0) ** 3)
        return Q, P, wss

    # solve at resting and stress (vasodilated) perfusion; WSS scales with flow, pressure drop too
    states = [("resting", PERF_REST, "#2c7fb8", "s--"), ("stress", PERF_STRESS, "#c0392b", "o-")]
    sol = {name: solve(q * mass_g / 60.0) for name, q, _, _ in states}
    diam = 2.0 * rad
    orders = np.asarray(orders)
    ords = sorted(int(o) for o in set(int(x) for x in orders) if o > 0)
    dm = np.array([np.median(diam[orders == o]) for o in ords])

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].axhspan(10, 70, color="0.88", label="physiological arterial (10-70)")   # Malek et al. 1999
    for name, q, col, mk in states:
        Q, P, wss = sol[name]
        wm = [np.median(wss[orders == o]) for o in ords]
        ax[0].plot(dm, wm, mk, color=col, lw=2, label=f"{name} ({q} mL/min/g)")
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("vessel diameter (mm)"); ax[0].set_ylabel("wall shear stress (dyn/cm2)")
    ax[0].set_title("Wall shear stress vs diameter"); ax[0].legend(fontsize=8)
    for name, q, col, mk in states:
        Q, P, wss = sol[name]
        pmed = [np.median(P[orders == o]) for o in ords]
        ax[1].plot(dm, pmed, mk, color=col, lw=2, label=f"{name} ({q} mL/min/g)")
    ax[1].set_xscale("log"); ax[1].set_xlabel("vessel diameter (mm)"); ax[1].set_ylabel("pressure (mmHg)")
    ax[1].set_title("Pressure vs diameter (down the tree)"); ax[1].legend(fontsize=8)
    fig.suptitle("Hemodynamics (lumped-parameter solve): inlet 100 mmHg, "
                 f"resting {PERF_REST} vs stress {PERF_STRESS} mL/min/g")
    fig.tight_layout(); fig.savefig(prefix + "_hemodynamics.png", dpi=150, bbox_inches="tight"); plt.close(fig)

    with open(prefix + "_hemodynamics.csv", "w", newline="") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["myocardial_mass_g", f"{mass_g:.2f}"])
        w.writerow(["resting_perfusion_mL/min/g", f"{PERF_REST}"])
        w.writerow(["stress_perfusion_mL/min/g", f"{PERF_STRESS}"])
        w.writerow(["inlet_pressure_mmHg", f"{Pin:.0f}"])
        w.writerow([])
        w.writerow(["order", "median_diameter_mm",
                    "WSS_resting_dyn/cm2", "WSS_stress_dyn/cm2",
                    "pressure_resting_mmHg", "pressure_stress_mmHg",
                    "flow_stress_mL/min"])
        Qr, Pr, wr = sol["resting"]; Qs, Ps, ws = sol["stress"]
        for k, o in enumerate(ords):
            idx = orders == o
            w.writerow([o, f"{dm[k]:.3f}",
                        f"{np.median(wr[idx]):.1f}", f"{np.median(ws[idx]):.1f}",
                        f"{np.median(Pr[idx]):.1f}", f"{np.median(Ps[idx]):.1f}",
                        f"{np.median(Qs[idx]) * 60.0:.4f}"])
    ws_all = sol["stress"][2]
    print(f"[hemo] lumped solve at resting {PERF_REST} and stress {PERF_STRESS} mL/min/g, "
          f"inlet 100 mmHg; stress per-order WSS "
          f"{min(np.median(ws_all[orders == o]) for o in ords):.0f}-"
          f"{max(np.median(ws_all[orders == o]) for o in ords):.0f} dyn/cm2 "
          f"-> {prefix}_hemodynamics.png/.csv")


def flow_calibrated_diameters(flow, orders, offset):
    """Option A (CCOV3): assign each segment a diameter as a single monotone function
    of its flow, calibrated so the per-order geometric-mean diameter matches the Kassab
    reference (MU). Because diameter stays an increasing function of flow, the tree is
    hemodynamically consistent (unlike the order-wise affine remap it replaces)."""
    flow = np.asarray(flow, float)
    orders = np.asarray(orders, int)
    xk, yk = [], []                                   # (log geo-mean flow, log Kassab diameter) per order
    for o in sorted(set(int(v) for v in orders if v > 0)):
        q = flow[orders == o]; q = q[q > 0]
        if q.size == 0:
            continue
        r = max(1, min(o + offset, len(MU)))
        xk.append(float(np.mean(np.log(q))))
        yk.append(float(np.log(MU[r - 1])))
    xk = np.asarray(xk); yk = np.asarray(yk)
    si = np.argsort(xk); xk = xk[si]; yk = yk[si]
    lq = np.log(np.maximum(flow, 1e-30))
    if xk.size < 2:
        ld = np.full_like(lq, yk[0] if xk.size else 0.0)
    else:
        ld = np.interp(lq, xk, yk)                     # monotone log-log interpolation
        s_lo = (yk[1] - yk[0]) / (xk[1] - xk[0])       # linear extrapolation beyond the knots
        s_hi = (yk[-1] - yk[-2]) / (xk[-1] - xk[-2])
        ld = np.where(lq < xk[0], yk[0] + s_lo * (lq - xk[0]), ld)
        ld = np.where(lq > xk[-1], yk[-1] + s_hi * (lq - xk[-1]), ld)
    return np.exp(ld)


def compute_literature_defined_orders(SO_init, diameters,
                                      max_iterations=100, convergence_threshold=0.005):
    """Diameter-defined Strahler ordering with the SAME update rule as the original method,
    but with the order boundaries fixed to the literature (Kassab MU/SD) instead of recomputed
    from the tree each iteration. It starts from the connection-based Strahler orders, then
    re-assigns every element by the literature thresholds using the original boundary formula
    (midpoint of the lower order's mean+SD and the upper order's mean-SD). Real literature
    order numbers are kept (no relabeling to 1..N), so a vessel is labeled by its true size."""
    MUk = np.asarray(MU, float); SDk = np.asarray(SD, float)
    DB = (MUk[1:] + SDk[1:] + MUk[:-1] - SDk[:-1]) / 2.0     # 10 boundaries between the 11 orders
    nB = MUk[-1] + SDk[-1]
    n_ord = len(MUk)
    di = np.asarray(diameters, float)
    SOdd = np.asarray(SO_init).astype(int)
    for _ in range(max_iterations):
        new = np.where(di > nB, n_ord + 1, n_ord - (di[:, None] < DB).sum(axis=1)).astype(int)
        change = float(np.mean(new != SOdd))
        SOdd = new
        if change < convergence_threshold:
            break
    return SOdd


def literature_order_of(diam):
    """Literature (Kassab) diameter-defined order for each diameter, using the original
    boundary formula on the reference means/SDs. Clipped to orders 1..11."""
    MUk = np.asarray(MU, float); SDk = np.asarray(SD, float)
    DB = (MUk[1:] + SDk[1:] + MUk[:-1] - SDk[:-1]) / 2.0
    o = len(MUk) - (np.asarray(diam, float)[:, None] < DB).sum(axis=1)
    return np.clip(o, 1, len(MUk)).astype(int)


def strahler_orders_per_segment(segment_connections, n):
    """Topological Strahler order per segment: leaves = 1; at a junction the parent takes the
    largest child order, +1 if two children share that largest order."""
    parent = segment_connections['parent']; children = segment_connections['children']
    strahler = np.zeros(n, int)
    roots = [i for i in range(n) if len(parent[i]) == 0]
    post = []; seen = set()
    for root in roots:
        st = [(root, False)]
        while st:
            nd, done = st.pop()
            if done: post.append(nd); continue
            if nd in seen: continue
            seen.add(nd); st.append((nd, True))
            for c in children[nd]: st.append((c, False))
    for s in post:
        ch = children[s]
        if not ch:
            strahler[s] = 1
        else:
            cor = [strahler[c] for c in ch]; mx = max(cor)
            strahler[s] = mx + 1 if cor.count(mx) >= 2 else mx
    return strahler


def group_strahler_elements(strahler, segment_connections, n):
    """Group segments into maximal connected runs of the SAME topological Strahler order
    (the proper 'element' definition). Returns seg2elem and the list of segment-id chains."""
    parent = segment_connections['parent']; children = segment_connections['children']
    seg2elem = np.full(n, -1, int); elem_segs = []
    for s in range(n):
        if seg2elem[s] != -1: continue
        o = strahler[s]; e = len(elem_segs); stk = [s]; seg2elem[s] = e; comp = [s]
        while stk:
            cur = stk.pop()
            for nb in parent[cur] + children[cur]:
                if seg2elem[nb] == -1 and strahler[nb] == o:
                    seg2elem[nb] = e; stk.append(nb); comp.append(nb)
        elem_segs.append(comp)
    return seg2elem, elem_segs


def refine_root_element(seg2elem, elem_segs, elem_order, diameters, segment_connections):
    """Robustness guard for the largest vessels. An element is labelled by the mean diameter of
    its segments; in a strongly asymmetric tree the top Strahler element can be a long trunk that
    tapers across a literature size boundary, so its mean dips below the order-11 cutoff and the
    root is mislabelled (e.g. order 10 instead of 11). If (and only if) the ROOT element spans
    more than one literature order, split it into maximal connected same-order runs so the biggest
    vessels always read their true order. No-op for a size-homogeneous root element (the usual
    symmetric case, e.g. gamma=3.0), leaving every other element untouched."""
    parent = segment_connections['parent']; children = segment_connections['children']
    n = len(seg2elem)
    roots = [i for i in range(n) if len(parent[i]) == 0]
    if not roots:
        return seg2elem, elem_segs, np.asarray(elem_order)
    e0 = int(seg2elem[roots[0]]); segs0 = list(elem_segs[e0])
    lit_of = {s: int(o) for s, o in zip(segs0, literature_order_of(diameters[segs0]))}
    if len(set(lit_of.values())) <= 1:
        return seg2elem, elem_segs, np.asarray(elem_order)           # homogeneous -> unchanged
    seg2elem = seg2elem.copy(); elem_segs = list(elem_segs); elem_order = list(map(int, elem_order))
    segset = set(segs0); visited = set(); first = True
    for s in segs0:
        if s in visited: continue
        o = lit_of[s]; stk = [s]; comp = []
        while stk:
            cur = stk.pop()
            if cur in visited: continue
            visited.add(cur); comp.append(cur)
            for nb in parent[cur] + children[cur]:
                if nb in segset and nb not in visited and lit_of[nb] == o: stk.append(nb)
        if first:
            elem_segs[e0] = comp; elem_order[e0] = o
            for x in comp: seg2elem[x] = e0
            first = False
        else:
            ne = len(elem_segs); elem_segs.append(comp); elem_order.append(o)
            for x in comp: seg2elem[x] = ne
    return seg2elem, elem_segs, np.asarray(elem_order)


def element_connectivity_counts(seg2elem, elem_order, segment_connections, n):
    """Directed element connectivity in the format expected by compute_mean_se:
    {(daughter_order m, parent_order n): [count per parent element of order n]}."""
    children = segment_connections['children']
    Ne = len(elem_order)
    e_children = [set() for _ in range(Ne)]
    for s in range(n):
        for c in children[s]:
            a, b = int(seg2elem[s]), int(seg2elem[c])
            if a != b:
                e_children[a].add(b)
    present = sorted(int(o) for o in set(int(x) for x in elem_order))
    counts = {(m, npar): [] for m in present for npar in present}
    for e in range(Ne):
        npar = int(elem_order[e])
        child_orders = [int(elem_order[c]) for c in e_children[e]]
        for m in present:
            counts[(m, npar)].append(child_orders.count(m))
    return counts


def plot_diameter_vs_order(prefix, orders, diameters):
    """The single figure produced per run: vessel diameter vs diameter-defined order, with
    our tree (median + interquartile band per order) overlaid on the literature (Kassab
    mean +/- SD). Log scale on diameter."""
    orders = np.asarray(orders); diameters = np.asarray(diameters, float)
    lit_ord = np.arange(1, len(MU) + 1)
    ords = sorted(int(o) for o in set(orders) if o > 0)
    med = np.array([np.median(diameters[orders == o]) for o in ords])
    lo  = np.array([np.percentile(diameters[orders == o], 25) for o in ords])
    hi  = np.array([np.percentile(diameters[orders == o], 75) for o in ords])

    fig = plt.figure(figsize=(8, 5.5))
    plt.plot(lit_ord, MU, 'o-', color='#2c3e50', lw=2, ms=6, label='Literature (Kassab)')
    plt.fill_between(lit_ord, MU - SD, MU + SD, color='#2c3e50', alpha=0.12)
    plt.plot(ords, med, 's--', color='#c0392b', lw=2, ms=7, label='Our tree (median)')
    plt.fill_between(ords, lo, hi, color='#c0392b', alpha=0.15, label='Our tree IQR')
    plt.yscale('log')
    plt.xlabel('Diameter-defined Strahler order'); plt.ylabel('Vessel diameter (mm)')
    plt.title('Diameter vs order: our tree vs literature')
    plt.xticks(range(1, len(MU) + 1)); plt.grid(True, which='both', ls=':', alpha=0.5); plt.legend()
    plt.tight_layout()
    plt.savefig(prefix + "_diameter_vs_order.png", dpi=150, bbox_inches='tight')
    plt.close(fig)


def run(xml_path, prefix, domain_vol=None):
    outdir = os.path.dirname(os.path.abspath(prefix)) or "."
    os.makedirs(outdir, exist_ok=True)

    # 1. Load geometry directly from the tree XML (no Swift / .dat / intermediate .vtp)
    vtk_polydata = build_polydata_from_tree_xml(xml_path)
    diameters = get_vessel_diameters(vtk_polydata)

    # 2-6. Proper two-step diameter-defined ordering with LITERATURE thresholds:
    #   (1) topological Strahler order per segment (the branching rule);
    #   (2) group segments into maximal same-Strahler-order elements;
    #   (3) relabel each element by its diameter against the fixed literature (Kassab)
    #       thresholds -- real order numbers kept, tree keeps its own diameters (no transform).
    segment_connections = build_segment_connections(vtk_polydata)
    ncell = vtk_polydata.GetNumberOfCells()
    strahler = strahler_orders_per_segment(segment_connections, ncell)
    seg2elem, elem_segs = group_strahler_elements(strahler, segment_connections, ncell)
    elem_order = literature_order_of([np.mean(diameters[c]) for c in elem_segs])
    # guard: keep the largest (root) vessels on their true order even if an asymmetric top trunk
    # averages below a size boundary; no-op for the usual symmetric (e.g. gamma=3.0) tree
    seg2elem, elem_segs, elem_order = refine_root_element(
        seg2elem, elem_segs, elem_order, diameters, segment_connections)
    SOF = np.array([elem_order[seg2elem[s]] for s in range(ncell)])
    max_tree_order = int(SOF.max())

    # 7. The single figure: diameter vs order, our tree overlaid on the literature
    plot_diameter_vs_order(prefix, SOF, diameters)

    # 8. Morphometry CSVs (connectivity computed on the SAME Strahler elements)
    connectivity_counts = element_connectivity_counts(seg2elem, elem_order, segment_connections, ncell)
    mean_se = compute_mean_se(connectivity_counts)
    unique_orders = np.unique(SOF)
    write_matrix_to_csv(mean_se, unique_orders, filename=prefix + "_connectivity_matrix.csv")
    write_diameter_stats_to_csv(compute_diameter_stats_per_order(SOF, diameters),
                                filename=prefix + "_diameter_stats.csv")
    write_diameter_ranges_to_csv(compute_diameter_ranges_per_order(SOF, diameters),
                                 filename=prefix + "_diameter_ranges.csv")

    # 9. Final vtp: transformed radius + final_order
    vtk_order = vtk.vtkIntArray(); vtk_order.SetName("final_order"); vtk_order.SetNumberOfComponents(1)
    vtk_order.SetNumberOfTuples(vtk_polydata.GetNumberOfCells())
    for i in range(vtk_polydata.GetNumberOfCells()):
        vtk_order.SetValue(i, int(SOF[i]))
    vtk_polydata.GetCellData().AddArray(vtk_order)
    vtk_polydata.GetCellData().SetScalars(vtk_order)
    w = vtk.vtkXMLPolyDataWriter(); w.SetFileName(prefix + "_tree.vtp")
    w.SetInputData(vtk_polydata); w.Write()

    # 10. Final obj: tube mesh of the transformed tree
    obj_poly = ensure_point_radius(vtk_polydata, "radius")
    tubes = build_tubes(obj_poly, "radius", n_sides=24, capping=True, fallback_radius=0.1)
    ow = vtk.vtkOBJWriter(); ow.SetFileName(prefix + "_tree.obj")
    ow.SetInputData(tubes); ow.Write()

    # 11. Hemodynamics: lumped-parameter resistive-network solve (paper Section 2.3.4)
    if domain_vol:
        compute_hemodynamics(xml_path, prefix, domain_vol, SOF)

    print(f"[postprocess] root order={max_tree_order}. Wrote: "
          f"{prefix}_tree.vtp, {prefix}_tree.obj, {prefix}_diameter_vs_order.png, "
          f"{prefix}_hemodynamics.png, {prefix}_connectivity_matrix.csv, "
          f"{prefix}_diameter_stats.csv, {prefix}_diameter_ranges.csv")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python3 coronary_postprocess.py tree.xml [output_prefix] [domain.vol]")
    xml = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(xml)[0]
    domain = sys.argv[3] if len(sys.argv) > 3 else None
    run(xml, prefix, domain)
