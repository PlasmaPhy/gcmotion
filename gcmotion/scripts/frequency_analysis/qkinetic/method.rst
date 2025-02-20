########################
# qkinetic calculation #
########################

For a profile with a specific mu, Pzeta and E, there could exist 0 or more
families of orbits, each represented by a contour line in the θ-Pθ phase
space. Each family corresponds to a differnet qkinetic.

We first check that the contour lines are valid, meaning the are completely
contained in the phase diagram and not get cutoff. Invalid contour lines are
discarded.

To calculate qkinetic for all the families, we iterate through all their
respective contour lines, and apply to each one the following steps:

1.  We create 2 new Profiles, with the same mu and E, but slightly lower and
slightly higher Pzetas (ideally infinitesimally close).

2.  We create the main ContourOrbit object, which is an object containing the
contour line's vertices, as well as some flags describing the orbit and the
methods to calculate them.

    2a. We calculate the Orbit's bounding box (e.g. the smallest rectangle
    fully containing the Orbit).

    2b. We create 2 new, "local" contours from each "adjacent" profile, larger
    than the bounding box, and try to find contour lines with the same E.

    2c. This again must yield at least 1 orbit in each new contour.  We must
    pick a correct size for the new contours to ensure that. We pick the
    closest ones by comparing the distances of their bounding boxes' bottom
    left corners. This method works for both trapped and passing orbits, and
    eliminates the errors deriving from the smallest orbits bounding boxes
    being comparable to the grid spacing.

2.  We create 2 new Profiles, with the same mu and E, but slightly lower and
slightly higher Pzetas (ideally infinitesimally close).

3.  For each, we create a new θ-Pθ contour, and we extract the contour lines
that correspond to the *same* energy.

4.  This will yield again, 0 or more contour lines. However, only one of them
corresponds to the main Orbit.

    3a. In the case that 0 lines are found, it means that the family no longer
    exists. This can either happen in the upper or lower pzeta contour, so we
    use the other one to calculate qkin. We also raise a flag that this is the
    final contour line we found.

    3b. In the case that 1 or more lines are found, we check again that the new
    lines are valid. If they are not, we end the calculation, since that means
    that the line was very close to be invalid in the first place.

4.  We create a ContourObrit object for each line, which apart from the line's
vertices, it also stores flags that classify the orbit as trapped,
co-/counter-passing, etc.

    4a. We first need to classify the orbit as trapped/passing, to correctly
    add the bottom points [(-2π, 0), (2π, 0)]

5.  We can easily calculate each line's bounding box (e.g. smallest rectangle
containing the whole line), and keep the closest one by comparing their bottom
left points.

    Note: Unfortunately, we can't do this for all the lines at the same time,
    since in the case that one family disappears, this method will choose the
    next closer contour line, which would be wrong. We can use the same grid,
    though.

6.  Now, for each line, we have 2 extra lines with slightly different Pzeta,
and therefore slightly different Jtheta. We classify them as trapped/passing in
order to add the 2 bottom points [(-2π, 0), (2π, 0)] and correctly calculate
their areas using the shoelace algorithm. Since the area is calculated on the
θ-Pθ phase space, it is by definition the corresponding Jθ of that orbit.

7.  We now have 2 very close contour obrits, with slightly different Pζ and
slightly different Jθ. We use those values to calculate q=dPζ/dJθ.

8.  We store qkin in the ContourObrit and return.

