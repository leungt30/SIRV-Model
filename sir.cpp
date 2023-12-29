#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "cpgplot.h"
using namespace std;

float birthrate; 
float deathrate;
float immunity_days;
double beta; 
double gamma_;
double vac_rate;
double vac_duration;
double vac_efficiency;

double dsdt(int population, double s, double i, double r, double v){
	// beta = disease transmission rate
	// - people infected + people born + people becoming vulnerable after immunity - people getting vaccinated + people whos vaccination expired
	return -beta*(s / population)*i + birthrate*population + r/immunity_days - vac_rate*s + v/vac_duration;
}
double didt(int population, double s, double i, double r, double v){
	// + infecting more people - people beating the disease - people dying from disease + vaccinated people infected
	return beta*(s/population)*i - gamma_*i - (deathrate * i) + beta*i*(1-vac_efficiency)*v/population;
} 
double drdt(int population, double s, double i, double r, double v){
	return gamma_*i - r/immunity_days;
}
double dvdt(int population, double s, double i, double r, double v){
	// vaccinated people - people whos vaccinations have worn off - people who got vac but still managed to get infected 
	return vac_rate * s - v/vac_duration - beta*i*(1-vac_efficiency)*v/population;
	// current model only says susceptible people get vaccination
}


int main(){
	double t_start = 0;
	double t_end;
	int n;
	double dt;
	int population;
	
	double days=1;
	int infected_start=0;
	// cout << "Enter number of steps, population, starting infected population, disease transmission rate, days to recover" << endl;
	// cin >> n >> population >> infected_start >> beta >> days ;
	ifstream inputfile;
	string filename = "config.txt"; //write all the inputs in here
	inputfile.open(filename.c_str());

	inputfile >> t_end >> n >> population >> infected_start >> beta >> days >> birthrate >> deathrate >> immunity_days >> vac_rate >> vac_duration >> vac_efficiency;

	//setup vars
	double t = t_start;
	double s = population-infected_start;
	double i = infected_start;
	double r = 0;
	double v = 0;

	double d = 0;
	gamma_ = 1.0 / days;
	double max_pop = population;
	

	//create arrays
	float tp[n+1],sp[n+1],ip[n+1],rp[n+1],popCount[n+1],vp[n+1];
	dt = (t_end - t_start) / n;
	for (int x = 0; x<n;x++){
		tp[x]=t;
		sp[x]=s;
		ip[x]=i;
		rp[x]=r;
		vp[x]=v;
		popCount[x]=s+i+r+v;
		if (max_pop < popCount[x]){
			max_pop = popCount[x];
		}
		
		double ks1 = dt * dsdt(popCount[x],s,i,r,v);
		double ks2 = dt * dsdt(popCount[x],s+0.5 * ks1, i,r,v);
		double ks3 = dt * dsdt(popCount[x],s+0.5 * ks2,i,r,v);
		double ks4 = dt * dsdt(popCount[x],s+ks3,i,r,v);
		
		double ki1 = dt * didt(popCount[x],s,i,r,v);
		double ki2 = dt * didt(popCount[x],s,i+0.5 * ki1,r,v);
		double ki3 = dt * didt(popCount[x],s,i+0.5 * ki2,r,v);
		double ki4 = dt * didt(popCount[x],s,i+ki3,r,v);
	
		double kr1 = dt * drdt(popCount[x],s,i,r,v);
		double kr2 = dt * drdt(popCount[x],s,i,r+0.5 * kr1,v);
		double kr3 = dt * drdt(popCount[x],s,i,r+0.5 * kr2,v);
		double kr4 = dt * drdt(popCount[x],s,i,r+kr3,v);
		
		double kv1 = dt * dvdt(popCount[x],s,i,r,v);
		double kv2 = dt * dvdt(popCount[x],s,i,r,v+0.5 * kv1);
		double kv3 = dt * dvdt(popCount[x],s,i,r,v+0.5 * kv2);
		double kv4 = dt * dvdt(popCount[x],s,i,r,v+kv3);

		t = t + dt;

		s = s+ (ks1+2*ks2+2*ks3+ks4)/6;
		i = i+ (ki1+2*ki2+2*ki3+ki4)/6;
		r = r+ (kr1+2*kr2+2*kr3+kr4)/6;
		v = v+ (kv1+2*kv2+2*kv3+kv4)/6;
		
		//cout << "total: " << s + i + r <<endl;

	}
	tp[n]=t;
	sp[n]=s;
	ip[n]=i;
	rp[n]=r;
	vp[n]=v;
	popCount[n]=s+i+r+v;

	// cout << dsdt(popCount[n],s,i,r,beta,gamma_) << endl << sp[n]-sp[n-5] << endl;
	cout << "R_0: " << beta/gamma_ << endl;
	




	//plotting
	if (!cpgopen("/XWINDOW")){ return 1;}
	
	//axis
	cpgenv(t_start,t_end,0.,max_pop*1.1,0,1);
	cpglab("t", "# of People", "SIRV Model");
	
	// susceptible
	cpgsci(3); //green 
	cpgline(n+1,tp,sp);
	// infected
	cpgsci(4); //blue estimate
	cpgline(n+1,tp,ip);
	// recovered
	cpgsci(5); //light blue 
	cpgline(n+1,tp,rp);
	// population count
	cpgsci(6); // purple
	cpgline(n+1,tp,popCount);//population count
	// vaccinated
	cpgsci(7); // yellow
	cpgline(n+1,tp,vp);//vaccination count


	cpgclos();
}
