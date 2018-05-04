#include <ros/ros.h>
#include <geometry_msgs/Twist.h>
#include <wiringPi.h>
#include <softPwm.h>
#include <math.h>


void callback(const geometry_msgs::Twist &msg)
{
	wiringPiSetupGpio();
	int pwm1=4, pwm2=26,direction1=14,direction2=6;
	softPwmCreate(pwm1, 0, 100);
	softPwmCreate(pwm2, 0, 100);
	
	pinMode(direction1, OUTPUT);
	pinMode(direction2, OUTPUT);
	if (msg.linear.x > 0.0)
	{
        std::cout << "hello" << std::endl;
		digitalWrite(direction1, HIGH);
		digitalWrite(direction2, LOW);		
	}
	else if (msg.linear.x < 0.0)
	{
		digitalWrite(direction1, LOW);
		digitalWrite(direction2, HIGH);
	}
	else if (msg.angular.z >= 0.0 && msg.linear.x == 0.0)
	{
		digitalWrite(direction1, LOW);
		digitalWrite(direction2, LOW);
	}
	else if (msg.angular.z < 0.0 && msg.linear.x == 0.0)
	{
		digitalWrite(direction1, HIGH);
		digitalWrite(direction2, HIGH);
	}
	float maxvel;
    std::cout << msg.linear.x << std::endl;
	maxvel = std::max(fabs(msg.linear.x), fabs(msg.angular.z));
    std::cout << maxvel << std::endl;

	softPwmWrite(pwm1, std::min(2.5f, maxvel) * 40.0f);
	softPwmWrite(pwm2, std::min(2.5f, maxvel) * 40.0f);
}
int main(int argc,char **argv)
{
	
	ros::init(argc, argv, "cmd_vel_subsriber");
	ros::NodeHandle n;
	ros::Subscriber sub = n.subscribe("/cmd_vel", 1000, callback);
	ros::Rate r(50);
	while (n.ok())
	{
		ros::spin();
	}
	return 0;
}
