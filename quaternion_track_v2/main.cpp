#include <stdio.h>
#include "csvReader.h"
#include "quaternion.h"

typedef struct sIMURec
{
	float t;
	float ax;
	float ay;
	float az;
	float wx;
	float wy;
	float wz;
}sIMURec;

class CIMUprocessor
{
private:
	sQ Qsg; //quaternion that translates current IMU orientation into ground frame,
	//with gravity vector equal to (0, 0, -1) and unknown x,y orientation. x,y orientation
	//can be kept too in case of IMU with magnetic sensor onboard
	sV cur_acceleration; //in ground frame
	sV cur_speed; //in ground frame
	sV cur_position; //in ground frame
public:
	CIMUprocessor()
	{
		Qsg.w = 1;
		Qsg.x = 0;
		Qsg.y = 0;
		Qsg.z = 0;
		cur_acceleration.x = 0;
		cur_acceleration.y = 0;
		cur_acceleration.z = 0;
		cur_speed.x = 0;
		cur_speed.y = 0;
		cur_speed.z = 0;
		cur_position.x = 0;
		cur_position.y = 0;
		cur_position.z = 0;
	}
	~CIMUprocessor();
	void update_g_vector(sV *acceleration, float coeff) //adjust quaternion in a way that given acceleration vector will mean (0, 0, -1) direction
	{ //coeff defines how much update will be there: 1 means that old position is ignored, and new one is set to given vector
	//0 means that new position is ignored, old one is kept
	//for regular updates values between 0.0001 and 0.01 make sense, experiments required to find optimal
		sV current_accel; //local variable in order not to change incoming value
		current_accel.x = acceleration->x;
		current_accel.y = acceleration->y;
		current_accel.z = acceleration->z;
		v_renorm(&current_accel); //normalize: we want only direction, not value

		sV target_accel; //we want acceleration to be translated into (0, 0, -1)
		target_accel.x = 0;
		target_accel.y = 0;
		target_accel.z = -1;
		
		//quaternion rotating V into U has form of [cos(theta/2), sin(theta/2)*(VxU)]
		//so we find it this way:
		float cos_theta = v_dot(&current_accel, &target_accel);
		float half_cos = sqrt(0.5f * (1.f + cos_theta));
		float half_sin = sqrt(0.5f * (1.f - cos_theta));
		sV update_axis = v_mult(&current_accel, &target_accel); //axis is cross product of vectors

		sQ update_q; //quaternion that defines updating rotation
		update_q.w = half_cos;
		update_q.x = half_sin*update_axis.x;
		update_q.y = half_sin*update_axis.y;
		update_q.z = half_sin*update_axis.z;
		q_renorm(&update_q); //normalize before use! :)

		//now we slowly adjust current quaternion into target one. For coeff less than 0.01,
		//this method gives good results. For higher adjustment speed, more complicated
		//method is required
		Qsg.w *= (1.0 - coeff);
		Qsg.x *= (1.0 - coeff);
		Qsg.y *= (1.0 - coeff);
		Qsg.z *= (1.0 - coeff);
		Qsg.w += coeff*update_q.w;
		Qsg.x += coeff*update_q.x;
		Qsg.y += coeff*update_q.y;
		Qsg.z += coeff*update_q.z;
		q_renorm(&Qsg); //normalize the result in order to keep numerical errors low
	};
	void update_rotation(sV *angular_speed_rads, double dt)
	{
		//angular speed means that our current quaternion should be rotated
		//in its own frame at angles equal to speed*dt
		sQ wQ; //angular speed quaternion

		wQ.w = 0;
		wQ.x = angular_speed_rads->x;
		wQ.y = angular_speed_rads->y;
		wQ.z = angular_speed_rads->z;
		
		sQ dQsg = q_mult(&Qsg, &wQ);
		Qsg.w += dQsg.w*0.5*dt;
		Qsg.x += dQsg.x*0.5*dt;
		Qsg.y += dQsg.y*0.5*dt;
		Qsg.z += dQsg.z*0.5*dt;
		q_renorm(&Qsg); //due to non-perfect calculations, resulting quaternion absolute value
						//might be slightly different from 1 - so we renormalize it
	};
	void update_position(sV *acceleration, double dt)
	{
		//simple integration in ground frame
		cur_acceleration.x = acceleration->x;
		cur_acceleration.y = acceleration->y;
		cur_acceleration.z = acceleration->z;
		rotate_v(&Qsg, &cur_acceleration); //rotate from sensor frame to ground frame
		cur_acceleration.z += 9.8; //compensate gravity
		
		cur_speed.x += dt*cur_acceleration.x;
		cur_speed.y += dt*cur_acceleration.y;
		cur_speed.z += dt*cur_acceleration.z;
		
		cur_position.x += dt*cur_speed.x;
		cur_position.y += dt*cur_speed.y;
		cur_position.z += dt*cur_speed.z;
	};
	void process_point(sV *acceleration, sV *angular_speed_rads, double dt)
	{
		double g_norm = v_norm(acceleration);
		if(g_norm > 0.9*9.8 && g_norm < 1.1*9.8)
			update_g_vector(acceleration, 0.001);
		
		update_rotation(angular_speed_rads, dt);
		update_position(acceleration, dt);
	};
	void rotate_to_ground(sV *v)
	{
		rotate_v(&Qsg, v);
	};
	sV get_position() {return cur_position;}
	sV get_speed() {return cur_speed;}
	sV get_acceleration() {return cur_acceleration;}
	void set_zero_position()
	{
		cur_acceleration.x = 0;
		cur_acceleration.y = 0;
		cur_acceleration.z = 0;
		cur_speed.x = 0;
		cur_speed.y = 0;
		cur_speed.z = 0;
		cur_position.x = 0;
		cur_position.y = 0;
		cur_position.z = 0;		
	};
};

int main(int argc, char **argv)
{
	csvReader *R;
	if(argc < 2)
	{
		printf("no input file specified - assuming ../accel_mag_o.csv \n");
		R = new csvReader("../accel_mag_o.csv");
	}
	else
		R = new csvReader(argv[1]);
	int cnt = R->getLinesCount();
	sIMURec *recs = new sIMURec[cnt];
	for(int x = 0; x < cnt; x++)
	{
		recs[x].t = R->readDouble();
		recs[x].ax = R->readDouble();
		recs[x].ay = R->readDouble();
		recs[x].az = R->readDouble();
		recs[x].wx = R->readDouble();
		recs[x].wy = R->readDouble();
		recs[x].wz = R->readDouble();
	}
	delete R;
	
	CIMUprocessor *imu_proc = new CIMUprocessor();
	
	float d2r = M_PI / 180.0;
	int g_timeout = 100;
	int g_low_time = g_timeout;
	int g_high_time = 0;
	
	sV cur_accel;
	sV cur_ang_speed;
	
	cur_accel.x = recs[0].ax * 9.81;
	cur_accel.y = recs[0].ay * 9.81;
	cur_accel.z = recs[0].az * 9.81;
	imu_proc->update_g_vector(&cur_accel, 1); //set initial position
	
	float event_max_speed = 0;
	int events_count;
	float max_speeds[1000]; //for now it's ok to assume that file has less than 1000 swings :)
	float max_times[1000];
	
	for(int x = 1; x < cnt-50; x++)
	{
		float dt = recs[x].t - recs[x-1].t;
		
		float gg = recs[x+30].ax*recs[x+30].ax + recs[x+30].ay*recs[x+30].ay + recs[x+30].az*recs[x+30].az;
		gg = sqrt(gg);

		if(gg >= 5.0)
		{
			g_high_time++;
			if(g_high_time > 5)
				g_low_time = 0;
		}
		
		if(gg < 5.0)
		{
			g_low_time++;
			g_high_time = 0;
			if(g_low_time > g_timeout)
			{
				imu_proc->set_zero_position();
				event_max_speed = 0;
			}
		}
		cur_accel.x = recs[x].ax * 9.81;
		cur_accel.y = recs[x].ay * 9.81;
		cur_accel.z = recs[x].az * 9.81;
		
		cur_ang_speed.x = d2r * recs[x].wx;
		cur_ang_speed.y = d2r * recs[x].wy;
		cur_ang_speed.z = d2r * recs[x].wz;
		
		imu_proc->process_point(&cur_accel, &cur_ang_speed, dt);
		
		if(1)if(g_low_time < g_timeout)
		{
			float vv = v_norm(imu_proc->get_speed());
			printf("%g, A %g %g %g %g ", recs[x].t, v_norm(imu_proc->get_acceleration()), imu_proc->get_acceleration().x, imu_proc->get_acceleration().y, imu_proc->get_acceleration().z);
			printf("V %g %g %g %g ", vv, imu_proc->get_speed().x, imu_proc->get_speed().y, imu_proc->get_speed().z);
			printf("P %g %g %g\n", imu_proc->get_position().x, imu_proc->get_position().y, imu_proc->get_position().z);
			if(vv > event_max_speed)
			{
				event_max_speed = vv;
				max_times[events_count] = recs[x].t;
			}
		}
		if(1)if(g_low_time == g_timeout)
		{
			max_speeds[events_count] = event_max_speed;
			events_count++;
			printf("\n\n");
		}

		if(0)if(x%20 == 1) //debug
		{
			{
				printf("%g, A %g %g %g %g\n", recs[x].t, v_norm(imu_proc->get_acceleration()), imu_proc->get_acceleration().x, imu_proc->get_acceleration().y, imu_proc->get_acceleration().z);
		/*		sV sx, sy, sz;
				sx.x = 1;
				sx.y = 0;
				sx.z = 0;
				sy.x = 0;
				sy.y = 1;
				sy.z = 0;
				sz.x = 0;
				sz.y = 0;
				sz.z = 1;
				rotate_v(&Qsg, &sx);
				rotate_v(&Qsg, &sy);
				rotate_v(&Qsg, &sz);

				printf("t: %g, S -> gnd:\n", recs[x].t); 
				printf("%g %g %g\n", sx.x, sx.y, sx.z);
				printf("%g %g %g\n", sy.x, sy.y, sy.z);
				printf("%g %g %g\n", sz.x, sz.y, sz.z);
				printf("accel: %g %g %g\n", accel.x, accel.y, accel.z);*/
			}
		}
	}
	
	float ms2mph = 2.237;
	for(int e = 0; e < events_count; e++)
	{
		printf("N %d time %g: max speed %g m/s (%g mph)\n", e, max_times[e], max_speeds[e], max_speeds[e]*ms2mph);
	}
	printf("mph values:\n");
	for(int e = 0; e < events_count; e++)
	{
		printf("%g\n", max_speeds[e]*ms2mph);		
	}
	
	return 0;
}
