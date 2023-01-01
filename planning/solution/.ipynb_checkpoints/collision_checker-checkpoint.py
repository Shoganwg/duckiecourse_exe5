import itertools
import random
from typing import List

from aido_schemas import Context, FriendlyPose
from dt_protocols import (
    Circle,
    CollisionCheckQuery,
    CollisionCheckResult,
    MapDefinition,
    PlacedPrimitive,
    Rectangle,
)

#################
import numpy as np
##################

__all__ = ["CollisionChecker"]


class CollisionChecker:
    params: MapDefinition

    def init(self, context: Context):
        context.info("init()")

    def on_received_set_params(self, context: Context, data: MapDefinition):
        context.info("initialized")
        self.params = data

    def on_received_query(self, context: Context, data: CollisionCheckQuery):
        collided = check_collision(
            environment=self.params.environment, robot_body=self.params.body, robot_pose=data.pose
        )
        result = CollisionCheckResult(collided)
        context.write("response", result)


def check_collision(environment: List[PlacedPrimitive], robot_body: List[PlacedPrimitive],
        robot_pose: FriendlyPose) -> bool:
    # This is just some code to get you started, but you don't have to follow it exactly

    # You can start by rototranslating the robot_body by the robot_pose
    rototranslated_robot: List[PlacedPrimitive] = []
    # == WRITE ME ==
    for ele in robot_body:
        ele.pose = robot_pose
        rototranslated_robot.append(ele)
    
    # Then, call check_collision_list to see if the robot collides with the environment
    collided = check_collision_list(rototranslated_robot, environment)
    
    return collided
    # # return a random choice
    # return random.uniform(0, 1) > 0.5


def check_collision_list(rototranslated_robot: List[PlacedPrimitive], environment: List[PlacedPrimitive]) -> bool:
    # This is just some code to get you started, but you don't have to follow it exactly
    for robot, envObject in itertools.product(rototranslated_robot, environment):
        if check_collision_shape(robot, envObject):
            return True

    return False

def get_dist_two_points(wcor1,wcor2):
    # wcor1 the world coordinate of point1 as a tuple, (float, float) of x and y value
    return ( (wcor1[0] - wcor2[0])**2  + (wcor1[1] - wcor2[1])**2   ) ** (1/2)

def is_two_circle_intersect(wc1x,wc1y,r1, wc2x,wcy2,r2):
    d12 = get_dist_two_points((wc1x,wc1y), (wc2x,wcy2))
    return d12 <= (r1+r2) 

def rec_to_world(rec:PlacedPrimitive):
    wrec = {}
    
    lps = [np.array([[rec.primitive.xmax], [rec.primitive.ymax]]),
           np.array([[rec.primitive.xmin], [rec.primitive.ymax]]),
           np.array([[rec.primitive.xmin], [rec.primitive.ymin]]),
           np.array([[rec.primitive.xmax], [rec.primitive.ymin]])]
    pose_at_w = np.array([[rec.pose.x],[rec.pose.y]])
    
    thera_r = np.radians(rec.pose.theta_deg)
    c, s = np.cos(thera_r), np.sin(thera_r)
    R = np.array(((c, -s), (s, c)))
    # wp = c + R@lp
    wvs = [(pose_at_w + R@lp).flatten() for lp in lps]
    wrec["wvs"] = wvs
    
    geo_c_at_w = ((wvs[0] + wvs[2]) / 2).flatten()
    wrec["center"] = geo_c_at_w
    wrec["rec_r1"] = min((rec.primitive.xmax - rec.primitive.xmin)/2, (rec.primitive.ymax - rec.primitive.ymin)/2)
    wrec["rec_r2"] = np.sqrt(((wvs[0] - wvs[2])**2).sum())/2
    
    return wrec

def cir_vertex(cirx, ciry, cir_r, rec_theta_deg):
    lp = np.array([[cir_r],[0]])
    # print(c1.x, c1.y)
    cir_cent_at_w = np.array([[cirx],[ciry]])
    Rs = []

    for deg_step in [0,90,180,270]:
        thera_r = np.radians(rec_theta_deg + deg_step)
        c, s = np.cos(thera_r), np.sin(thera_r)
        R = np.array(((c, -s), (s, c)))
        Rs.append(R)
    # wp = c + R@lp
    wvs = [(cir_cent_at_w + Rn@lp).flatten() for Rn in Rs]
    return wvs 

def is_same_side_of_line(lA,lB,pC,pD):
    AB = lB - lA
    AC = pC - lA
    AD = pD - lA
    return np.dot(np.cross(AB,AC), np.cross(AB,AD)) >= 0

def is_point_in_rec(rec1vs, rec1c, point):
    rec1vs.append(rec1vs[0])
    in_count = 0
    for i in range(4):
        lA,lB,pC,pD = rec1vs[i],rec1vs[i+1],rec1c,point
        if is_same_side_of_line(lA,lB,pC,pD):
            in_count += 1
            if in_count == 4:
                return True

    return False
        
    
def is_any_rec_vertex_in_circle(rec1vs, wc2x,wcy2,r2):
    for vertex in rec1vs:
        dist = get_dist_two_points((vertex[0],vertex[1]), (wc2x,wcy2))
        if dist <= r2:
            return True

    return False

    
def is_two_rec_intersect(rec1vs, rec1c, rec2vs, rec2c):
 
    for i in range(4):
        if is_point_in_rec(rec2vs, rec2c, rec1vs[i]):
            return True
        if is_point_in_rec(rec1vs, rec1c, rec2vs[i]):
            return True
        
    
    return False
    
def check_collision_shape(a: PlacedPrimitive, b: PlacedPrimitive) -> bool:
    # This is just some code to get you started, but you don't have to follow it exactly

    # This is just some code to get you started, but you don't have to follow it exactly
    if isinstance(a.primitive, Circle) and isinstance(b.primitive, Circle):
        
        # == WRITE ME ==      
        return is_two_circle_intersect(a.pose.x, a.pose.y, a.primitive.radius , b.pose.x, b.pose.y, b.primitive.radius)
        
    elif isinstance(a.primitive, Rectangle) and isinstance(b.primitive, Circle):
        # == WRITE ME ==
        wrec = rec_to_world(a)
        cvts = cir_vertex(b.pose.x, b.pose.y, b.primitive.radius, a.pose.theta_deg)
        if not is_two_circle_intersect(wrec["center"][0], wrec["center"][1], wrec["rec_r2"] , b.pose.x, b.pose.y, b.primitive.radius):
            return False
        
        if is_any_rec_vertex_in_circle(wrec["wvs"], b.pose.x, b.pose.y, b.primitive.radius):
            return True
        
        if is_two_rec_intersect(wrec["wvs"], wrec["center"], cvts, np.array([b.pose.x, b.pose.y])):
            return True


    elif isinstance(a.primitive, Circle) and isinstance(b.primitive, Rectangle):
        # == WRITE ME ==
        wrec = rec_to_world(b)
        cvts = cir_vertex(a.pose.x, a.pose.y, a.primitive.radius, b.pose.theta_deg)
        if not is_two_circle_intersect(wrec["center"][0], wrec["center"][1], wrec["rec_r2"] , a.pose.x, a.pose.y, a.primitive.radius):
            return False
        
        else:
            return is_two_rec_intersect(wrec["wvs"], wrec["center"], cvts, np.array([a.pose.x, a.pose.y]))


        
    elif isinstance(a.primitive, Rectangle) and isinstance(b.primitive, Rectangle):
        # == WRITE ME ==
        wreca = rec_to_world(a)
        wrecb = rec_to_world(b)
        return is_two_rec_intersect(wreca["wvs"], wreca["center"], wrecb["wvs"], wrecb["center"])
        
    # ...
    # # for now let's return a random guess

    else: 
        raise Exception(f"primitive combination of this pair never been considere before. {a.primitive} {b.primitive}") 
