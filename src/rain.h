/**
 * Rain Simulation Module
 * 
 * This header simulates rain by feeding a liquid film within a specified region.
 * The impact region is defined by a line with a Gaussian distribution perpendicular to it,
 * creating a "rain band" effect.
 * 
 * Line equation: normal.x*x + normal.y*y = intercept
 * 
 * The Gaussian distribution is oriented along the line (tangent direction),
 * allowing rain intensity to vary smoothly across the impact zone.
 */

attribute {
    struct rain {
        coord u;          // Incoming velocity of rain
        double rvf;       // Rain volume fraction

        double stddev;    // Standard deviation of Gaussian impact distribution (width of region)
        double mean;      // Mean position of Gaussian along the line (in tangent direction)

        coord normal;     // Unit normal vector to the line (points to from rain region)
        coord tangent;    // Unit tangent vector along the line
        double intercept; // Line intercept parameter (α in n·x = α)

        coord pref;       // Reference point on the line for coordinate system origin
    } rain;
}

extern scalar f;

/**
 * Initialize reference point and tangent vector for the rain line.
 * 
 * The reference point is set to a convenient intercept:
 *   - x-intercept (α/n_x, 0) if n_x ≠ 0
 *   - y-intercept (0, α/n_y) if n_y ≠ 0
 * 
 * The tangent vector is perpendicular to the normal, oriented along the line.
 */
event defaults (i = 0)
{
    if (f.rain.normal.x != 0.)
        f.rain.pref = (coord) {f.rain.intercept / f.rain.normal.x, 0}; // use x-intercept
    else if (f.rain.normal.y != 0.) 
        f.rain.pref = (coord) {0, f.rain.intercept / f.rain.normal.y}; // use y-intercept

    // Tangent is perpendicular to normal: rotate normal by 90° counterclockwise
    f.rain.tangent = (coord){-f.rain.normal.y, f.rain.normal.x};
}

/**
 * Evaluate unnormalized Gaussian distribution.
 * 
 * Returns exp(-0.5 * ((h - mean) / stddev)²)
 * 
 * @note The probability = 100% at the top of the bell curve (since its unnormalized)
 * 
 * @param stddev Standard deviation of the distribution
 * @param mean   Center position of the distribution
 * @param h      Position coordinate along the line
 * @return       Gaussian weight at position h
 */
double gauss_distribution(double stddev, double mean, double h)
{
    return exp(-0.5*sq((h - mean)/stddev));
}

/** 
 * Check if a point is inside the rain impact region.
 * 
 * Uses the signed distance to the line to determine which side the point is on.
 * The normal vector points away from the rain region, so:
 *   - Positive distance: inside rain region (gets rain)
 *   - Negative distance: outside rain region (no rain)
 * 
 * @param c  Scalar field containing line parameters
 * @param pc Point coordinates to check
 * @return   true if point is in rain region, false otherwise
 */
bool is_in_impact_region(scalar c, coord pc)
{
    // Compute signed distance: n·p - α
    return c.rain.normal.x*pc.x + c.rain.normal.y*pc.y - c.rain.intercept >= 0;
}

/**
 * Calculate rain probability at a given point.
 * 
 * Algorithm:
 *   1. Project the point onto the rain line
 *   2. Compute the coordinate h along the line (tangent direction) from reference point
 *   3. Evaluate Gaussian distribution at h
 *   4. Return 0 if point is outside the impact region
 * 
 * @param point Current grid point (unused but required by interface)
 * @param pc    Point coordinates in space
 * @param c     Scalar field containing rain parameters
 * @return      Rain probability weight (0 to 1)
 */
double rain_probability(Point point, coord pc, scalar c)
{
    // Step 1: Project point onto the line
    // Compute signed distance t from point to line
    double tt = c.rain.normal.x*pc.x + c.rain.normal.y*pc.y - c.rain.intercept;

    // Project by moving distance t in the normal direction: p_proj = p - t*n
    coord pc_proj = {pc.x - tt*c.rain.normal.x, pc.y - tt*c.rain.normal.y};

    // Step 2: Compute coordinate along the line
    // h is the signed distance from reference point along tangent direction
    double h = (pc_proj.x - c.rain.pref.x)*(c.rain.tangent.x) + (pc_proj.y - c.rain.pref.y)*(c.rain.tangent.y);

    // Step 3: Return Gaussian weight if inside region, 0 otherwise
    return is_in_impact_region(c, pc)? gauss_distribution(c.rain.stddev, c.rain.mean, h): 0;
}
