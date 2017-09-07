import numpy
import MDAnalysis
import spin
print "here is spin"
print dir(spin)

def spinangle_from_trajectory(psf, trajectory, spinref, selection):
    """ 
    Extracts the spinangle of the selection over the trajectory.
    Uses MDAnalysis and a homemade C module. Could be upgraded.
    psf is the name of the topology PSF file, trajectory the DCD file, spinref the PDB that has the reference strcture and selection 
    is a string with a MDAnalysis selection.
    MDAnalysis can also replace the PSF with a PDB and the DCD with a XTC. Untested with Gromacs trajectories.
    
    Example arguments:
    psf = "/media/ub/DATA/Miscellaneous/init_AqpZ_prots.psf"
    trajectory = "/media/ub/DATA/Trajectories/ABFa/ABFa0-143_prots.unwrap.dcd"
    spinref = "/media/ub/DATA/Miscellaneous/spinDimer_reference.pdb"
    selection = "name BAS and bynum 475:948"
    """
    
    # Load the trajectory
    mobile = MDAnalysis.Universe(psf, trajectory)
    # Get the coordinates from this trajectory. This step could be upgraded.
    coords_nocenter = extractCoordsMDa(mobile, selection)
    # Load the reference and get the coordinates
    reference = MDAnalysis.Universe(spinref)
    reference_nocenter = reference.select_atoms(selection).positions
    
    # Center the coordinates on the protein - we want the spin angle, just remove the centroid
    # Get contiguous arrays using numpy, necessary for the C function
    coords_center = []
    for i in range(len(coords_nocenter)):
        coords_center.append(numpy.ascontiguousarray(shift_center(coords_nocenter[i]), dtype=numpy.float64))
    reference_center = numpy.ascontiguousarray(shift_center(reference_nocenter), dtype=numpy.float64)

    # Call get_spinangle_traj
    spin_timeserie = []
    for frame in coords_center:
        spinangle = spin.get_spinangle_traj(frame, reference_center, len(reference_center))
        spin_timeserie.append(spinangle)
    # Return a numpy array
    return numpy.array(spin_timeserie)

def shift_center(conformation):
    """Center and typecheck the conformation"""

    conformation = numpy.asarray(conformation)
    if not conformation.ndim == 2:
        raise ValueError('conformation must be two dimensional')
    _, three = conformation.shape
    if not three == 3:
        raise ValueError('conformation second dimension must be 3')

    centroid = numpy.mean(conformation, axis=0)
    centered = conformation - centroid
    return centered

def extractCoordsMDa(mobile, selection):
    """ Extracts the coordinates with MDanalysis """

    # Use MDAnalysis built-in C routines to get the timeserie
    coords = MDAnalysis.collection.addTimeseries(MDAnalysis.Timeseries.Atom('v', mobile.select_atoms(selection)))
    MDAnalysis.collection.compute(mobile.trajectory)
    tmp = []
    # This is embarassing - the collection object has all the coordinates in the wrong order, so we have to use a dirty thing to get the right one
    # I should write a C thing for this, or we should be able to do something using numpy
    for frame in range(len(MDAnalysis.collection[0][0][0])):
        tmp.append([[MDAnalysis.collection[0][i][0][frame], MDAnalysis.collection[0][i][1][frame], 
        MDAnalysis.collection[0][i][2][frame]] for i in range(len(MDAnalysis.collection[0]))])
        if frame == 5:
            break

    # Yay
    MDAnalysis.collection.clear()

    # Give back a numpy array
    return numpy.array(tmp)
