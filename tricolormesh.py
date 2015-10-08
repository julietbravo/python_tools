import numpy as np
from matplotlib.collections import PolyCollection
import pylab as plt

def tricolormesh(x_vert, y_vert, z=None, vmin=None, vmax=None, edgecolor='none', linewidth=0.2, cmap=None, rasterized=False):
    """
    Plot (filled) triangular mesh
    Author: Bart van Stratum, Max-Planck Institute for Meteorology, Oct 2015

    Keyword arguments:
    x_vert     -- 2 dimensional np.array {cells, edges} with x location of cell edges
    y_vert     -- 2 dimensional np.array {cells, edges} with y location of cell edges
    z          -- np.array of size(cells) with values per cell
    vmin       -- minimal value of z 
    vmax       -- maximum value of z 
    edgecolor  -- color of the cell edges (or 'none' to disable edges)
    linewidth  -- linewidth of cell edges
    cmap       -- colormap for cell colors
    rasterized -- switch to enable/disable rasterization 
    """
    
    # Set default color map
    if cmap is None:
        cmap = plt.cm.jet

    # Some (very limited..) input checks
    if(x_vert.ndim != 2 or y_vert.ndim != 2):
        raise Exception('x_vert and y_vert must have 2 dimensions [cells,edges_of_cell]')

    # Clip the z array to take vmin and vmax into account
    if z is not None:
        z_clip = z.copy()
        if(vmin is not None):
            z_clip[z<vmin] = vmin
        if(vmax is not None):
            z_clip[z>vmax] = vmax
    else:
        z_clip = z
    
    # Create iterator with correct format for input in PolyCollection()
    #tri = (zip(x_vert[i,:],y_vert[i,:]) for i in range(x_vert[:,0].size))
    tri = [zip(x_vert[i,:],y_vert[i,:]) for i in range(x_vert[:,0].size)]

    # Create the polygons
    if z is None:
        col = PolyCollection(tri, cmap=cmap, edgecolors=edgecolor, facecolors='none', linewidths=linewidth, rasterized=rasterized)
    else:
        col = PolyCollection(tri, array=z_clip, cmap=cmap, edgecolors=edgecolor, linewidths=linewidth, rasterized=rasterized)

    # Add to axes
    ax = plt.gca()  
    ax.add_collection(col)

    # Draw second time such that there are really no edges visible (bit hackish...)
    if(edgecolor == 'none'):
        ax.add_collection(col)

    # Scale axis to wrap PolyCollection
    ax.autoscale_view() 

    # Return the polycollection, for drawing a colorbar
    return col


if __name__ == "__main__":
    import netCDF4 as nc4

    nc = nc4.Dataset('example_data/ICON_example_data.nc', 'r')
    clon  = nc.variables['clon' ][:]                                       
    clat  = nc.variables['clat' ][:]                                       
    z_ifc = np.random.random(clon.size)                                        
    clonv = nc.variables['clon_vertices'][:,:]                                        
    clatv = nc.variables['clat_vertices'][:,:]                                        
    nc.close()

    plt.figure()
    tricolormesh(clonv, clatv, z_ifc, edgecolor='k', linewidth=0.2, cmap=plt.cm.terrain, rasterized=False)
