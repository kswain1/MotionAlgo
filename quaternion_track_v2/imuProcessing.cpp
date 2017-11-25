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

