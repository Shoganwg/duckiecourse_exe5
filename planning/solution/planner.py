########################
import numpy as np
np.typeDict = np.sctypeDict
np.int = np.intc

from dt_protocols import (
    Circle,
    CollisionCheckQuery,
    CollisionCheckResult,
    MapDefinition,
    PlacedPrimitive,
    Rectangle,
)

from collision_checker import *
from collision_checker import check_collision_list, rec_to_world

import networkx as nx
########################


from typing import List

from aido_schemas import Context, FriendlyPose
from dt_protocols import (
    PlacedPrimitive,
    PlanningQuery, PlanningResult, PlanningSetup, PlanStep, Rectangle,
)

__all__ = ["Planner"]





class Planner:
    params: PlanningSetup

    def init(self, context: Context):
        context.info("init()")

    def on_received_set_params(self, context: Context, data: PlanningSetup):
        context.info("initialized")
        self.params = data

        # This is the interval of allowed linear velocity
        # Note that min_velocity_x_m_s and max_velocity_x_m_s might be different.
        # Note that min_velocity_x_m_s may be 0 in advanced exercises (cannot go backward)
        max_velocity_x_m_s: float = self.params.max_linear_velocity_m_s
        min_velocity_x_m_s: float = self.params.min_linear_velocity_m_s
        
        ##############################
        # max_angular_velocity_deg_s: float = self.params.max_angular_velocity_deg_s
        #############################
        
        # This is the max curvature. In earlier exercises, this is +inf: you can turn in place.
        # In advanced exercises, this is less than infinity: you cannot turn in place.
        max_curvature: float = self.params.max_curvature

        # these have the same meaning as the collision exercises
        body: List[PlacedPrimitive] = self.params.body
        environment: List[PlacedPrimitive] = self.params.environment

        # these are the final tolerances - the precision at which you need to arrive at the goal
        tolerance_theta_deg: float = self.params.tolerance_theta_deg
        tolerance_xy_m: float = self.params.tolerance_xy_m

        # For convenience, this is the rectangle that contains all the available environment,
        # so you don't need to compute it
        bounds: Rectangle = self.params.bounds
    
    def opt_turn(self, dg1, dg2):
        diffd = dg2 - dg1
        deltg= 0 # initialize degree to turn
        if diffd > 180:
            deltg = -(360-diffd)
        elif diffd < -180:
            deltg = 360-diffd
        else:
            deltg = diffd
        return deltg
    
    def gen_turn_step(self, turndg):
        # this only can be use when self.params.
        dur = abs(turndg)/self.params.max_angular_velocity_deg_s
        vel = 0.0
        angv = np.sign(turndg) * self.params.max_angular_velocity_deg_s
        
        return PlanStep(dur,vel,angv)
    
    def gen_straight_step(self, x0,y0, xt, yt):
        
        dist = ((xt-x0)** 2 + (yt-y0)**2) ** (1/2)

        dur = float(dist/self.params.max_linear_velocity_m_s)
        vel = float(self.params.max_linear_velocity_m_s)
        angv = 0.0
        
        return PlanStep(dur,vel,angv)
    
    def turn_straight_turn(self, FP1, FP2):
        x1,y1, thd1 ,x2, y2, thd2 = FP1.x, FP1.y, FP1.theta_deg ,FP2.x, FP2.y, FP2.theta_deg
        dirx,diry = x2-x1, y2-y1
        str_thd = np.degrees(np.arctan2(diry, dirx))

        turn1dg = self.opt_turn(thd1, str_thd)
        turn2dg = self.opt_turn(str_thd, thd2)
        
        plan = [self.gen_turn_step(turn1dg), self.gen_straight_step(x1,y1,x2, y2), self.gen_turn_step(turn2dg)]
        
        return plan
    
    def get_third_length(self, a, b, rad):
        return (a**2 + b**2 -2*a*b*np.cos(rad))**(1/2)
    
    def if_circle_in_zone(self,PlPr):
        
        r = PlPr.primitive.radius
        
        xmi,xma,ymi,yma = PlPr.pose.x-r,PlPr.pose.x+r,PlPr.pose.y-r,PlPr.pose.y+r
        too_high = yma>= self.params.bounds.ymax
        too_low  = ymi<= self.params.bounds.ymin
        too_right= xma>= self.params.bounds.xmax
        too_left = xmi<= self.params.bounds.ymin
        
        if too_high or too_low or too_right or too_left:
            return False
        else:
            return True
    
    def get_min_circle_zone_radius(self):
        # rbw = rec_to_world(rb_plpr)
        rb_plpr = self.params.body[0]
        cvr = 1 / self.params.max_curvature # min curvature radius 
        
        xma,yma,xmi,ymi = rb_plpr.primitive.xmax,rb_plpr.primitive.ymax,rb_plpr.primitive.xmin,rb_plpr.primitive.ymin
        px,py = rb_plpr.pose.x, rb_plpr.pose.y
        
        l1 = (xma**2 + yma**2)**(1/2)
        l2 = (xmi**2 + yma**2)**(1/2)
        l3 = (xmi**2 + ymi**2)**(1/2)
        l4 = (xma**2 + ymi**2)**(1/2)
        
        rad1 = np.pi/2 + np.arctan2(xma-px,yma-py)
        rad2 = (3/2)*np.pi - np.arctan2(xmi-px,yma-py)
        rad3 = np.pi + np.arctan2(xmi-px,ymi-py) + np.pi/2
        rad4 = np.pi - np.arctan2(xma-px,ymi-py)
        
        r1 = self.get_third_length(l1,cvr,rad1)
        r2 = self.get_third_length(l2,cvr,rad2)
        r3 = self.get_third_length(l3,cvr,rad3)
        r4 = self.get_third_length(l4,cvr,rad4)
                
        r = max(r1,r2,r3,r4)   # it robot rotating, what is the radius of the min circle 
        
        return r
    
    def get_circle_from_plpr(self, rb_plpr, clock_deg=1.0):
            
        pose_at_w = np.array([[rb_plpr.pose.x],[rb_plpr.pose.y]])
        r_cuvr = 1 / self.params.max_curvature # min curvature radius
        c_local_ctcw = np.array([[0],[r_cuvr]])
        c_local_cw = np.array([[0],[-r_cuvr]])
        theta_deg = rb_plpr.pose.theta_deg

        thera_r = np.radians(theta_deg)
        c, s = np.cos(thera_r), np.sin(thera_r)
        R = np.array(((c, -s), (s, c)))
        # wp = c + R@lp
        c_ctcw_x,c_ctcw_y = (pose_at_w + R@c_local_ctcw).flatten()   # if robot run counter clockwise, what is the center of the circle
        c_cw_x,c_cw_y = (pose_at_w + R@c_local_cw).flatten()
        
        r = self.get_min_circle_zone_radius()
        
        # theta_deg = 1.0 means it is counter clockwise, if 2.0 means it is clockwise
        pp_ctcw = PlacedPrimitive(pose=FriendlyPose(x=c_ctcw_x,
                                                    y=c_ctcw_y,
                                                    theta_deg=clock_deg), 
                                    primitive=Circle(radius=r))
        
        # pp_cw = PlacedPrimitive(pose=FriendlyPose(x=c_cw_x,
        #                                             y=c_cw_y,
        #                                             theta_deg=clock_deg), 
        #                             primitive=Circle(radius=r))
        
        print("clock_deg == 1.0", clock_deg == 1.0)
        print("self.if_circle_in_zone(pp_ctcw)", self.if_circle_in_zone(pp_ctcw))
        print("(not check_collision_list([pp_ctcw], self.params.environment))",(not check_collision_list([pp_ctcw], self.params.environment)))
        
        # if clock_deg == 1.0 and self.if_circle_in_zone(pp_ctcw) and (not check_collision_list([pp_ctcw], self.params.environment)):
        if clock_deg == 1.0 and self.if_circle_in_zone(pp_ctcw):
            
            return True, pp_ctcw
        
        # elif clock_deg == 2.0 and self.if_circle_in_zone(pp_cw) and not check_collision_list([pp_cw], self.params.environment):
        #     return True,pp_cw
        
        else:
            # context.info(f"ATTENTION ############################ could not generate start_or_end circle ")
            print(f"ATTENTION ############################ could not generate start_or_end circle ")
            return False, None
        
        # if self.if_circle_in_zone(pp_cw) and not check_collision_list([pp_cw], self.params.environment):
        #     key2 = f"{start_or_end}_cw"
        #     circle_dict[key2] = pp_cw
            
        # if len(circle_dict) == 0:
        #     context.info(f"ATTENTION ############################ could not generate {start_or_end} circle ")
        

    def is_tract_doable(self, pp_tract_rec_list):   
        chk = check_collision_list(pp_tract_rec_list, self.params.environment)
        # print(f"The track result is ######### {chk}")
        # print(pp_tract_rec_list)
        # print(self.params.environment)
        return not check_collision_list(pp_tract_rec_list, self.params.environment)
    
    def get_track_recpp(self,distance,br_deg,rec_px, rec_py, tengantial_list):
        xma,xmi,yma,ymi = self.params.body[0].primitive.xmax, self.params.body[0].primitive.xmin,\
                        self.params.body[0].primitive.ymax,self.params.body[0].primitive.ymin
        
        x_width = self.params.body[0].primitive.xmax - self.params.body[0].primitive.xmin
        px1,py1,px2,py2 = tengantial_list
        
        
        ###############################
        nk = int(distance/(x_width))
        if nk<= 1:
            onepp = rec_pp = PlacedPrimitive(pose=FriendlyPose(x=(px1+px2)/2,
                                                        y=(py1+py2)/2,
                                                        theta_deg=br_deg), 
                                        primitive=Rectangle(xmax=xma,xmin=xmi,ymax=yma,ymin=ymi))
            return [onepp]
        else:
            rec_list = []
            for i in range(nk):
                rec_px = px1 + (px2-px1)*(i/(nk-1))
                rec_py = py1 + (py2-py1)*(i/(nk-1))
                rec_pp = PlacedPrimitive(pose=FriendlyPose(x=rec_px,
                                                            y=rec_py,
                                                            theta_deg=br_deg), 
                                            primitive=Rectangle(xmax=xma,xmin=xmi,ymax=yma,ymin=ymi))
                rec_list.append(rec_pp)

            return rec_list
    
        
    def bridge_two_circles(self, cpp1, cpp2):
        """
        return is_birdge_able, how long, connection_plan
        """
        rb_width = self.params.body[0].primitive.ymax - self.params.body[0].primitive.ymin
        deg_shift = {}
        deg_shift["counterclockwise"] = -90
        deg_shift["clockwise"] = 90
        
        x1,y1,x2,y2 = cpp1.pose.x, cpp1.pose.y, cpp2.pose.x, cpp2.pose.y
        
        br_deg = np.degrees(np.arctan2(y2-y1, x2-x1)) 
        if br_deg < 0.0:
            br_deg = br_deg + 360
                
        
        distance = ((x2-x1)**2 + (y2-y1)**2)**(1/2)
        dt = distance / self.params.max_linear_velocity_m_s
        step_plan = PlanStep(duration=dt, \
                         angular_velocity_deg_s=0.0,\
                         velocity_x_m_s=self.params.max_linear_velocity_m_s)
        
        rcur = 1.0/self.params.max_curvature
               
        if cpp1.pose.theta_deg == 1.0 and cpp1.pose.theta_deg == cpp2.pose.theta_deg :
            # theta_deg == 1, means this is a counter clockwise circle
            dx = rcur*np.cos(np.radians(br_deg + deg_shift["counterclockwise"]))
            dy = rcur*np.sin(np.radians(br_deg + deg_shift["counterclockwise"]))
            px1,py1,px2,py2 = x1+dx,y1+dy,x2+dx,y2+dy
            rec_px, rec_py = (px1+px2)/2, (py1+py2)/2
            rec_pp_list = self.get_track_recpp(distance,br_deg,rec_px, rec_py, [px1,py1,px2,py2])
                             
            if self.is_tract_doable(rec_pp_list):
                br_plan = PlanStep(duration=dt, 
                                   angular_velocity_deg_s=0.0,
                                   velocity_x_m_s=self.params.max_linear_velocity_m_s)
                return True, dt, br_deg, br_plan, [x1,y1,x2,y2]            
            else:
                return False, None, None, None, [x1,y1,x2,y2]
         
        else:
            return False, None, None, None
#         elif cpp1.pose.theta_deg == 2.0 and cpp1.pose.theta_deg == cpp2.pose.theta_deg :
#             dx = rcur*np.cos(np.radians(br_deg + deg_shift["clockwise"]))
#             dy = rcur*np.sin(np.radians(br_deg + deg_shift["clockwise"]))
#             px1,py1,px2,py2 = x1+dx,y1+dy,x2+dx,y2+dy
#             rec_px, rec_py = (px1+px2)/2, (py1+py2)/2

            

            
        # else:
        #     return False, None, None, None 
        
    def turn_within_circle_ctcw(self,br_deg1,br_deg2):
        vel = self.params.max_linear_velocity_m_s
        aglv = self.params.max_angular_velocity_deg_s
        
        if br_deg2 >= br_deg1:
            deltdg = br_deg2 - br_deg1
        elif br_deg2 < br_deg1:
            deltdg = 360 - (br_deg1 - br_deg2)    
        
        dt = deltdg/aglv
            
        turn_plan = PlanStep(duration=dt, 
                           angular_velocity_deg_s=aglv,
                           velocity_x_m_s=vel)
        return dt, turn_plan
    
    def generate_circle_from_xy(self,x,y, clock_deg):
        r = self.get_min_circle_zone_radius()
        pp = PlacedPrimitive(pose=FriendlyPose(x=x,
                                                y=y,
                                                theta_deg=clock_deg), 
                                    primitive=Circle(radius=r))
        if self.if_circle_in_zone(pp) and not check_collision_list([pp], self.params.environment):
             return True, pp
        else:
            return False, None
    
    def generate_valid_grid_circle(self, crl_plpr_0, clock_deg):
        """the most strick valid circle, 
        if motion, freeze and combime all obsturcles"""
        
        r = crl_plpr_0.primitive.radius
        
        bounds_pp = PlacedPrimitive(pose=FriendlyPose(x=0, y=0, theta_deg=0), primitive=self.params.bounds)
        
        bounds_wrec = rec_to_world(bounds_pp) 
        xmi, xma, ymi, yma = bounds_wrec["wvs"][2][0] + r, bounds_wrec["wvs"][0][0] -r, bounds_wrec["wvs"][2][1]+r, bounds_wrec["wvs"][0][1] -r
        
        away = 6*r
        
        xn, yn = int((xma-xmi)/away),int((yma-ymi)/away)
        
        existable_circle_dict = {}
        for i in range(xn+1):
            for j in range(yn+1):
                cx,cy = xmi+i*(2*r),ymi+j*(2*r)
                is_good, pp = self.generate_circle_from_xy(cx,cy, clock_deg)
                if is_good:
                    existable_circle_dict[(i,j)] = pp
        
        return existable_circle_dict
        # TO DO
        # When it is a motion cases
        
    
    def generate_best_plan(self, PQ, clock_deg=1.0):
        """return bool, plan"""
        p0 = PQ.start
        pt = PQ.target
        
        rb_plpr_0 = PlacedPrimitive(pose=p0,primitive=self.params.body[0].primitive)
        rb_plpr_t = PlacedPrimitive(pose=pt,primitive=self.params.body[0].primitive)
        
        good0, crl_plpr_0 = self.get_circle_from_plpr(rb_plpr_0, clock_deg)
        goodt, crl_plpr_t = self.get_circle_from_plpr(rb_plpr_t, clock_deg)
        
        if (not good0) or (not goodt):
            print(good0, crl_plpr_0)
            print(goodt, crl_plpr_t)
        #     context.info("##########################  can not generate  cir0 or cirt")
        
        transit_circles_dict = self.generate_valid_grid_circle(crl_plpr_0, clock_deg)
        
        rot = 180.0 / self.params.max_angular_velocity_deg_s
        
        G = nx.DiGraph()
        
        can0t, dt0t, br_deg0t, br_plan0t, x1y1x2y2_0t = self.bridge_two_circles(crl_plpr_0, crl_plpr_t) 
        if can0t:
            G.add_edge("0", "t", weight=rot+dt0t, planstep = br_plan0t, br_deg = br_deg0t, cord = x1y1x2y2_0t )
            
            deg0 = p0.theta_deg
            degt = pt.theta_deg
            
            path = nx.shortest_path(G, "0", "t")
            plan = []
            for i in range(len(path)):
                    # ['0', (0, 1), 't']
                    if (i == 0) and (path[i] == "0"):
                        d = G.get_edge_data('0', path[i+1])
                        # rotate
                        dt, turn_plan0 = self.turn_within_circle_ctcw(deg0, d["br_deg"]) 
                        plan.append(turn_plan0)
                        # straight
                        print(d["cord"])
                        plan.append(d["planstep"])

                    elif (i == (len(path) -1)) and (path[i] == "t"):

                        # rotate
                        d = G.get_edge_data(path[i-1], "t")
                        dt, turn_plant = self.turn_within_circle_ctcw(d["br_deg"], degt) 
                        print(d["cord"])
                        plan.append(turn_plant)
            print("Get it , return now ################################################################3")
            # context.info("Get it , return now !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            return True, plan
        
        elif not can0t:
            for trs_node, pp in transit_circles_dict.items():
                can01, dt01, br_deg01, br_plan01, x1y1x2y2_01 = self.bridge_two_circles(crl_plpr_0, pp) 
                can1t, dt1t, br_deg1t, br_plan1t, x1y1x2y2_1t = self.bridge_two_circles(crl_plpr_t, pp) 

                if can01:
                    G.add_edge("0", trs_node, weight=rot+dt01, planstep = br_plan01, br_deg = br_deg01, cord= x1y1x2y2_01 )

                if can1t:
                    G.add_edge(trs_node, "t", weight=rot+dt1t, planstep = br_plan1t, br_deg = br_deg1t, cord = x1y1x2y2_1t )


            
            trs_key_list = list(transit_circles_dict)
            kn = len(trs_key_list)
            for i in range(kn):
                for j in range(kn):
                    if j==i:
                        pass
                    if j!=i:
                        can, dt,br_deg, br_plan,x1y1x2y2 = self.bridge_two_circles(transit_circles_dict[trs_key_list[i]],\
                                                                                   transit_circles_dict[trs_key_list[j]])
                        if can:
                            G.add_edge(trs_key_list[i], trs_key_list[j], weight=rot+dt, planstep = br_plan, br_deg = br_deg, cord=x1y1x2y2)

            feasible = False

            path = []

            try:
                path = nx.shortest_path(G, "0", "t")
                feasible = True
            except:
                feasible = False
                return False, None
                
            print(path)

            if len(path) > 0:
                plan = []
                d01 = G.get_edge_data('0', path[1])
                deg0 = p0.theta_deg
                degt = pt.theta_deg

                for i in range(len(path)):
                        # ['0', (0, 1), 't']
                        if (i == 0) and (path[i] == "0"):
                            d = G.get_edge_data('0', path[i+1])
                            # rotate
                            dt, turn_plan0 = self.turn_within_circle_ctcw(deg0, d["br_deg"]) 
                            plan.append(turn_plan0)
                            # straight
                            print(d["cord"])
                            plan.append(d["planstep"])

                        elif (i == (len(path) -1)) and (path[i] == "t"):

                            # rotate
                            d = G.get_edge_data(path[i-1], "t")
                            dt, turn_plant = self.turn_within_circle_ctcw(d["br_deg"], degt) 
                            print(d["cord"])
                            plan.append(turn_plant)

                            # no straight
                        else:
                            dlast = G.get_edge_data(path[i-1], path[i])
                            dnext = G.get_edge_data(path[i], path[i+1])
                            # rotate
                            dt, turn_plan = self.turn_within_circle_ctcw(dlast["br_deg"], dnext["br_deg"]) 
                            plan.append(turn_plan)
                            # straight
                            print(dnext["cord"])
                            plan.append(dnext["planstep"])

                return feasible, plan
            
            else:
                return False, None
        return False,None    
            
    def on_received_query(self, context: Context, data: PlanningQuery):
        # A planning query is a pair of initial and goal poses
        start: FriendlyPose = data.start
        goal: FriendlyPose = data.target
        
        
        # You start at the start pose. You must reach the goal with a tolerance given by
        # tolerance_xy_m and tolerance_theta_deg.

        # You need to declare if it is feasible or not
        feasible = False

        if (len(self.params.environment) == 0) and (self.params.max_curvature >1000):
            plan = self.turn_straight_turn(start, goal)
            result: PlanningResult = PlanningResult(True, plan)
            print(f"#################feasible 1, return with plan len()= {len(plan)},\n {plan}")
            context.write("response", result)
            return result
        
        else:
            feasible, plan = self.generate_best_plan( data, clock_deg=1.0)
            if feasible:
                print(f"####################### feasible 2, return with plan len()= {len(plan)},\n {plan}")
                result: PlanningResult = PlanningResult(True, plan)
                context.write("response", result)
                return result
            
        
            if not feasible:
                # If it's not feasible, just return this.
                print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n infeasible, return with []  instead of None")
                result: PlanningResult = PlanningResult(False, [])
                context.write("response", result)
                return result
        
        
        
        # If it is feasible you need to provide a plan.

        # A plan is a list of PlanStep
        # plan: List[PlanStep] = []
        
#         class PlanningSetup(MapDefinition):
#         bounds: Rectangle
#         max_linear_velocity_m_s: float
#         min_linear_velocity_m_s: float
#         max_angular_velocity_deg_s: float
#         max_curvature: float
#         tolerance_xy_m: float
#         tolerance_theta_deg: float

 # environment: List[PlacedPrimitive]
 #    body: List[PlacedPrimitive]
        
        # if (len(self.params.environment) == 0) and (self.params.max_curvature >1000):
        #     plan = self.turn_straight_turn(start, goal)
        #     result: PlanningResult = PlanningResult(True, plan)
        #     context.info("response", result)
        # A plan step consists in a duration, a linear and angular velocity.
        
        # if feasible:
        #     # If it's not feasible, just return this.
        #     result: PlanningResult = PlanningResult(True, plan)
        #     context.write("response", result)
        #     return result
        
#         # Let's trace a square of side L at maximum velocity.
#         L = 1.0
#         duration_straight_m_s = L / self.params.max_linear_velocity_m_s
#         duration_turn_deg_s = 90.0 / self.params.max_angular_velocity_deg_s
#         # The plan will be: straight, turn, straight, turn, straight, turn, straight, turn

#         straight = PlanStep(duration=duration_straight_m_s, angular_velocity_deg_s=0.0,
#                             velocity_x_m_s=self.params.max_linear_velocity_m_s)
#         turn = PlanStep(duration=duration_turn_deg_s,
#                         angular_velocity_deg_s=self.params.max_angular_velocity_deg_s,
#                         velocity_x_m_s=0.0)

#         plan.append(straight)
#         plan.append(turn)
#         plan.append(straight)
#         plan.append(turn)
#         plan.append(straight)
#         plan.append(turn)
#         plan.append(straight)
#         plan.append(turn)

#         result: PlanningResult = PlanningResult(feasible, plan)
#         context.write("response", result)
