"""
bhtree compatibility module using scipy.spatial.cKDTree

This module provides a drop-in replacement for the legacy MGLTools bhtree
library using scipy's cKDTree for efficient neighbor searches.

The original bhtree is a compiled C extension that is difficult to install
on modern systems. This pure Python replacement provides the same API
using scipy's optimized KDTree implementation.
"""

import numpy as np
from scipy.spatial import cKDTree


class bhtreelib:
    """Compatibility wrapper for bhtree.bhtreelib"""
    
    class BHtree:
        """
        Ball-tree / BH-tree implementation using scipy's cKDTree.
        
        Provides efficient O(n log n) neighbor searches instead of O(n²)
        brute force, which is critical for large molecules like proteins.
        """
        
        def __init__(self, coords, radii=None, granularity=10):
            """
            Initialize the BHtree with atom coordinates and radii.
            
            Args:
                coords: List or array of (x, y, z) coordinates
                radii: List of atomic radii (bond order radii)
                granularity: Ignored (kept for API compatibility)
            """
            self.coords = np.asarray(coords, dtype=np.float64)
            self.radii = np.asarray(radii, dtype=np.float64) if radii is not None else None
            self.tree = cKDTree(self.coords)
        
        def closePointsPairsInTree(self, factor=1.175):
            """
            Find all pairs of atoms that are close enough to be bonded.
            
            Uses the sum of bond order radii * factor as the distance cutoff.
            
            Args:
                factor: Multiplier for sum of radii (default 1.175)
                
            Returns:
                List of (i, j) index pairs for atoms close enough to bond
            """
            if self.radii is None:
                # Use a default cutoff if no radii provided
                max_bond_dist = 2.0
            else:
                # Maximum possible bond distance
                max_radius = np.max(self.radii)
                max_bond_dist = 2 * max_radius * factor
            
            # Query all pairs within max possible bond distance
            pairs = self.tree.query_pairs(r=max_bond_dist, output_type='ndarray')
            
            if len(pairs) == 0:
                return []
            
            if self.radii is not None:
                # Filter pairs based on actual sum of radii
                valid_pairs = []
                for i, j in pairs:
                    actual_dist = np.linalg.norm(self.coords[i] - self.coords[j])
                    cutoff = (self.radii[i] + self.radii[j]) * factor
                    if actual_dist <= cutoff:
                        valid_pairs.append((i, j))
                return valid_pairs
            else:
                return [tuple(p) for p in pairs]
        
        def closePointsDist2(self, point, cutoff, indices_out, dist2_out):
            """
            Find all points within cutoff distance of a query point.
            
            Args:
                point: Query point (x, y, z)
                cutoff: Maximum distance
                indices_out: Output array for indices (modified in place)
                dist2_out: Output array for squared distances (modified in place)
                
            Returns:
                Number of neighbors found
            """
            point = np.asarray(point)
            
            # Query neighbors within cutoff
            indices = self.tree.query_ball_point(point, r=cutoff)
            n = len(indices)
            
            if n > 0:
                # Fill output arrays
                for i, idx in enumerate(indices[:len(indices_out)]):
                    indices_out[i] = idx
                    diff = self.coords[idx] - point
                    dist2_out[i] = np.dot(diff, diff)
            
            return min(n, len(indices_out))


# For direct import compatibility
BHtree = bhtreelib.BHtree
