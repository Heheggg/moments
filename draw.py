from display import *
from matrix import *
from math import *
from gmath import *

def add_polygon( polygons, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point(polygons, x0, y0, z0);
    add_point(polygons, x1, y1, z1);
    add_point(polygons, x2, y2, z2);

def draw_polygons( matrix, screen,zb, color ):
    if len(matrix) < 2:
        print 'Need at least 3 points to draw'
        return

    point = 0
    r = 255
    g = 255
    b = 255
    
    while point < len(matrix) - 2:

        normal = calculate_normal(matrix, point)[:]
        print normal
        if normal[2] > 0:
            scanline(matrix,screen,[r,g,b])
            draw_line( int(matrix[point][0]),
                       int(matrix[point][1]),
                       int(matrix[point+1][0]),
                       int(matrix[point+1][1]),
                       screen, color)
            draw_line( int(matrix[point+2][0]),
                       int(matrix[point+2][1]),
                       int(matrix[point+1][0]),
                       int(matrix[point+1][1]),
                       screen, color)
            draw_line( int(matrix[point][0]),
                       int(matrix[point][1]),
                       int(matrix[point+2][0]),
                       int(matrix[point+2][1]),
                       screen, color)    
            point+= 3

            r = (r + 31) % 256
            g = (r + 109) % 256
            b = (r + 199) % 256

        point += 3
        
def scanline(matrix, screen, color):
    sorty = {matrix[0][0]:0,matrix[0][1]:1,matrix[0][2]:2}
    output = sorted(sorty.keys())
    p_low = matrix[sorty[output[0]]]
    p_mid = matrix[sorty[output[1]]]
    p_top = matrix[sorty[output[2]]]
    

    y0 = p_low[1]
    x0 = p_low[0]
    d_x0 = ((1.0*(p_top[0]-p_low[0]))/(p_top[1]-p_low[1]))
    d_z0 = ((1.0*(p_top[2]-p_low[2]))/(1.0*(p_top[1]-p_low[1])))

    if p_low[1] != p_mid[1]:
        d_x1 = ((1.0*(p_mid[0]-p_low[0]))/(1.0*(p_mid[1]-p_low[1])))
        d_z1 = ((1.0*(p_mid[2]-p_low[2]))/(1.0*(p_mid[1]-p_low[1])))
        x1 = p_low[0]
        z1 = p_low[2]
    else:
        d_x1 = ((1.0*(p_top[0]-p_low[0]))/(1.0*(p_top[1]-p_low[1])))
        d_z1 = ((1.0*(p_top[2]-p_low[2]))/(1.0*(p_top[1]-p_low[1])))
        x1 = p_mid[0]
        z1 = p_mid[2]

    y = p_low[1]
    i = 0
    while y < p_top[1]:
        drawline(int(x0),int(y),int(z0),int(x1),int(y), int(z1), screen, zbuffer, color)

        if(y < p_mid[1] and p_mid[1]-y <1):
            
            x0 = p_low[0] + (p_mid[1]-p_low[1])*d_x0
            x1 = p_mid[0]
            y = p_mid[1]
            z0 = p_low[2] + (p_mid[1]-p_low[1])*d_z0
            z1 = p_mid[2]
            drawline(int(x0),int(y),int(z0),int(x1),int(y), int(z1), screen, zbuffer, color)

        if y == p_mid[1]:
            x1 = p_mid[0]
            z1 = p_mid[2]
            if p_top[1] != p_mid[1]:
                d_x1 = ((1.0*(p_top[0]-p_mid[0]))/(1.0*(p_top[1]-p_mid[1])))
                d_z1 = ((1.0*(p_top[2]-p_low[2]))/(1.0*(p_top[1]-p_mid[1])))
            else:
                d_x1 = ((1.0*(p_top[0]-p_low[0]))/(1.0*(p_top[1]-p_low[1])))
                d_z1 = ((1.0*(p_top[2]-p_low[2]))/(1.0*(p_top[1]-p_low[1])))

        x0 += d_x0
        x1 += d_x1
        y += 1
        z0 += d_z0
        z1 += d_z1

def add_box( polygons, x, y, z, width, height, depth ):
    x1 = x + width
    y1 = y - height
    z1 = z - depth

    #front
    add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
    add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);
  
    #back
    add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
    add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);
  
    #right side
    add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
    add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
    #left side
    add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
    add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);
  
    #top
    add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
    add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
    #bottom
    add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
    add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);

def add_sphere( edges, cx, cy, cz, r, step ):
    points = generate_sphere(cx, cy, cz, r, step)
    num_steps = int(1/step+0.1)
    
    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps

    num_steps+= 1
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):
            
            p0 = lat * (num_steps) + longt
            p1 = p0+1
            p2 = (p1+num_steps) % (num_steps * (num_steps-1))
            p3 = (p0+num_steps) % (num_steps * (num_steps-1))

            if longt != num_steps - 2:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p1][0],
		             points[p1][1],
		             points[p1][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2])
            if longt != 0:
	        add_polygon( edges, points[p0][0],
		             points[p0][1],
		             points[p0][2],
		             points[p2][0],
		             points[p2][1],
		             points[p2][2],
		             points[p3][0],
		             points[p3][1],
		             points[p3][2])

def generate_sphere( cx, cy, cz, r, step ):
    points = []
    num_steps = int(1/step+0.1)
    
    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps
            
    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop+1):
            circ = step * circle

            x = r * math.cos(math.pi * circ) + cx
            y = r * math.sin(math.pi * circ) * math.cos(2*math.pi * rot) + cy
            z = r * math.sin(math.pi * circ) * math.sin(2*math.pi * rot) + cz

            points.append([x, y, z])
            #print 'rotation: %d\tcircle%d'%(rotation, circle)
    return points
        
def add_torus( edges, cx, cy, cz, r0, r1, step ):
    points = generate_torus(cx, cy, cz, r0, r1, step)
    num_steps = int(1/step+0.1)
    
    lat_start = 0
    lat_stop = num_steps
    longt_start = 0
    longt_stop = num_steps
    
    for lat in range(lat_start, lat_stop):
        for longt in range(longt_start, longt_stop):

            p0 = lat * (num_steps) + longt;
            if (longt == num_steps - 1):
	        p1 = p0 - longt;
            else:
	        p1 = p0 + 1;
            p2 = (p1 + num_steps) % (num_steps * num_steps);
            p3 = (p0 + num_steps) % (num_steps * num_steps);

            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p3][0],
                        points[p3][1],
                        points[p3][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2] )
            add_polygon(edges,
                        points[p0][0],
                        points[p0][1],
                        points[p0][2],
                        points[p2][0],
                        points[p2][1],
                        points[p2][2],
                        points[p1][0],
                        points[p1][1],
                        points[p1][2] )

def generate_torus( cx, cy, cz, r0, r1, step ):
    points = []
    num_steps = int(1/step+0.1)
    
    rot_start = 0
    rot_stop = num_steps
    circ_start = 0
    circ_stop = num_steps

    print num_steps
    
    for rotation in range(rot_start, rot_stop):
        rot = step * rotation
        for circle in range(circ_start, circ_stop):
            circ = step * circle

            x = math.cos(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cx;
            y = r0 * math.sin(2*math.pi * circ) + cy;
            z = -1*math.sin(2*math.pi * rot) * (r0 * math.cos(2*math.pi * circ) + r1) + cz;

            points.append([x, y, z])
    return points

def add_circle( points, cx, cy, cz, r, step ):
    x0 = r + cx
    y0 = cy
    t = step

    while t <= 1.00001:
        x1 = r * math.cos(2*math.pi * t) + cx;
        y1 = r * math.sin(2*math.pi * t) + cy;

        add_edge(points, x0, y0, cz, x1, y1, cz)
        x0 = x1
        y0 = y1
        t+= step

def add_curve( points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type ):

    xcoefs = generate_curve_coefs(x0, x1, x2, x3, curve_type)[0]
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, curve_type)[0]

    t = step
    while t <= 1.00001:
        x = xcoefs[0] * t*t*t + xcoefs[1] * t*t + xcoefs[2] * t + xcoefs[3]
        y = ycoefs[0] * t*t*t + ycoefs[1] * t*t + ycoefs[2] * t + ycoefs[3]
                
        add_edge(points, x0, y0, 0, x, y, 0)
        x0 = x
        y0 = y
        t+= step

def draw_lines( matrix, screen, color ):
    if len(matrix) < 2:
        print 'Need at least 2 points to draw'
        return
    
    point = 0
    while point < len(matrix) - 1:
        draw_line( int(matrix[point][0]),
                   int(matrix[point][1]),
                   int(matrix[point+1][0]),
                   int(matrix[point+1][1]),
                   screen, color)    
        point+= 2
        
def add_edge( matrix, x0, y0, z0, x1, y1, z1 ):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)
    
def add_point( matrix, x, y, z=0 ):
    matrix.append( [x, y, z, 1] )
    
