import random

from display import *
from matrix import *
from math import *
from gmath import *

def magnitude(v):
    return math.sqrt(sum(i**2 for i in v))

def normalize(v):
    return [i/magnitude(v) for i in v]

def dot_product(a,b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def get_color(normal, lighting_info, lighting_name):
    ambient = get_ambient(lighting_info, lighting_name)
    diffuse = get_diffuse(normal, lighting_info, lighting_name)
    specular = get_specular(normal, lighting_info, lighting_name)
    I = [ambient[0] + diffuse[0] + specular[0],
         ambient[1] + diffuse[1] + specular[1],
        ambient[2] + diffuse[2] + specular[2]]
    return [int(max(min(x, 255), 0)) for x in I]

def add_polygon( polygons, x0, y0, z0, x1, y1, z1, x2, y2, z2 ):
    add_point(polygons, x0, y0, z0);
    add_point(polygons, x1, y1, z1);
    add_point(polygons, x2, y2, z2);

def draw_polygons( matrix, screen,zb, lighting_info, lighting_name  ):
    if len(matrix) < 2:
        print 'Need at least 3 points to draw'
        return

    point = 0
    while point < len(matrix) - 2:

        normal = calculate_normal(matrix, point)[:]
        print normal
        if normal[2] > 0:
            color = get_color(normal, lighting_info, lighting_name)
            scanline(matrix,point,screen,zb,color)
            
            draw_line( int(matrix[point][0]),
                       int(matrix[point][1]),
                       matrix[point][2],
                       int(matrix[point+1][0]),
                       int(matrix[point+1][1]),
                       matrix[point+1][2],
                       screen, zb, color)
            draw_line( int(matrix[point+2][0]),
                       int(matrix[point+2][1]),
                       matrix[point+2][2],
                       int(matrix[point+1][0]),
                       int(matrix[point+1][1]),
                       matrix[point+1][2],
                       screen, zb, color)
            draw_line( int(matrix[point][0]),
                       int(matrix[point][1]),
                       matrix[point][2],
                       int(matrix[point+2][0]),
                       int(matrix[point+2][1]),
                       matrix[point+2][2],
                       screen, zb, color)   

        point += 3
        
def scanline(matrix, i, screen, zbuffer, color):
    output = sorted([matrix[i],matrix[i+1],matrix[i+2]],key=lambda x:(x[1],x[0]))
    p_low = output[0]
    p_mid = output[1]
    p_top = output[2]
    

    z0 = p_low[2]
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
        draw_line(int(x0),int(y),int(z0),int(x1),int(y), int(z1), screen, zbuffer, color)

        if(y < p_mid[1] and p_mid[1]-y <1):
            
            x0 = p_low[0] + (p_mid[1]-p_low[1])*d_x0
            x1 = p_mid[0]
            y = p_mid[1]
            z0 = p_low[2] + (p_mid[1]-p_low[1])*d_z0
            z1 = p_mid[2]
            draw_line(int(x0),int(y),int(z0),int(x1),int(y), int(z1), screen, zbuffer, color)

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

def get_ambient(lighting_info, lighting_name):
    source_color = lighting_info["ambient"]
    const = lighting_info["constants"][lighting_name]
    ambient_const = [ const["red"][0], const["green"][0], const["blue"][0] ]
    ambient = [ ambient_const[0] * source_color[0], #red
                ambient_const[1] * source_color[1], #green
                ambient_const[2] * source_color[2] ] #blue
    return ambient
    
def get_diffuse(normal, lighting_info, lighting_name):
    const = lighting_info["constants"][lighting_name]
    diffuse_const = [ const["red"][1], const["green"][1], const["blue"][1] ]

    normal = normalize(normal)

    total_diffuse = [0, 0, 0]
    for entry in lighting_info["lights"]:
        light_source = lighting_info["lights"][entry]
        source_color = light_source["color"]
        location = light_source["location"]
        dot_prod = dot_product(normal, location)

        diffuse = [ source_color[0] * diffuse_const[0] * dot_prod,
                    source_color[1] * diffuse_const[1] * dot_prod,
                    source_color[2] * diffuse_const[2] * dot_prod ]
        total_diffuse = [ total_diffuse[0] + diffuse[0],
                          total_diffuse[1] + diffuse[1],
                          total_diffuse[2] + diffuse[2] ]
    return total_diffuse

def get_specular(normal, lighting_info, lighting_name):
    const = lighting_info["constants"][lighting_name]
    specular_const = [ const["red"][1], const["green"][1], const["blue"][1] ]

    normal = normalize(normal)

    total_specular = [0, 0, 0]
    for entry in lighting_info["lights"]:
        light_source = lighting_info["lights"][entry]
        source_color = light_source["color"]
        location = light_source["location"]
        dot_prod = dot_product(normal, location)

        a = [2 * dot_prod * x for x in normal]
        b = [ a[0] - location[0],
              a[1] - location[1],
              a[2] - location[2] ]
        b = normalize(b)
        view = normalize([1, 1, 1])
        c = dot_product(b, view)
        
        specular_r = source_color[0] * specular_const[0] * c
        specular_g = source_color[1] * specular_const[1] * c
        specular_b = source_color[2] * specular_const[2] * c

        total_specular = [ total_specular[0] + specular_r,
                           total_specular[1] + specular_g,
                           total_specular[2] + specular_b ]

    return total_specular

        
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
    
def draw_line( x0, y0, z0, x1, y1, z1, screen, zbuffer, color ):
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt

    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)
    wide = False
    tall = False

    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B

    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y)
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y

    while ( loop_start < loop_end ):
        plot( screen, zbuffer, color, x, y, z )
        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):
            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
        else:
            x+= dx_east
            y+= dy_east
            d+= d_east
        loop_start+= 1

    plot( screen, zbuffer, color, x, y, z )

