var atmos = {
	// Singleton class for transforming Standard Atmosphere variables
	// Author: Jimmy Lu
	
	Density: function(sigma, unit){
		// Calculate the density given density ratio
		//
		// Inputs:	sigma: density ratio
		//			unit: 0 for imperial (slug/cu.ft), 1 for SI (kg/cu.m)
		// Outputs:	density

		var rho0 = 0.00237689; // slug/cu.ft
		var density = rho0*sigma;

		if (unit==1){
			density = density*515.379;
		}
		
		return density;
	},
	
	PressureAltitude: function(delta, unit){
		// Calculate the pressure altitude given pressure ratio
		// Input: delta: pressure ratio
		//        unit: 0 = imperial (ft), 1 = SI (m)
		// Output: pressure altitude in ft or m

		var pressure_altitude = 0;
		if (delta > 0.2234)
			pressure_altitude = (1-Math.pow(delta,1/5.2559))/6.87535e-6;
		else
			pressure_altitude = 36089 - 20806*Math.log(4.477*delta);

		if (unit==1)
			pressure_altitude = pressure_altitude*0.3048;
		
		return pressure_altitude;
	},
	
	StaticPressure: function(delta, unit){
		// Calculate static pressure (ambient pressure) given delta
		//
		// Inputs: delta: pressure ratio
		//         unit: 0=imperial (psf), 1=SI (N/sq.m)
		// Outputs: static pressure in lbf/sq.ft or N/sq.m


		var P0 = 2116.224;
		var P = P0*delta;

		if (unit==1)
			P = P*47.88026;
		
		return P;
	},
	
	TotalPressure: function(P, qc){
		// Calculate total pressure given static pressure and impact pressure
		//
		// Inputs:	P: static pressure
		//			qc: impact pressure
		// Outputs:	total pressure

		var Pt = P + qc;
		return Pt;
	},
	
	SpeedOfSound: function(theta, unit){
		// Calculate the speed of sound given the temperature ratio
		// Inputs: theta: temp ratio
		//         unit: 0=imperial (kt), 1=SI (m/s)
		// Outputs: speed of sound in kt or m/s

		var a0 = 661.48;    // [kt]
		var a = a0*Math.pow(theta,0.5);

		if (unit==1)
			a = a*0.514444;

		return a;
	},
	
	OAT: function(theta){
		// Calculate the delta ISA and outside air temperature in C given
		// temperature ratio and pressure ratio
		//
		// Inputs: theta: temp ratio
		
		var OAT0 = 288.15; // Kelvin
		var OAT = theta*OAT0;
		
		// Convert to C
		OAT = OAT - 273.15;
		
		return OAT;
	},
	
	DISA: function(theta, delta){
		// Calculate the delta ISA and outside air temperature in C given
		// temperature ratio and pressure ratio
		//
		// Inputs: theta: temp ratio
		//         delta: pressure ratio

		// Convert to K
		var C2K = 273.15;
		var OAT = this.OAT(theta) + C2K;
		var OAT0 = this.OAT(1) + C2K;
		
		// ISA temperature at delta:
		var theta_ISA = this.Theta(delta, 0);
		var DISA = OAT - OAT0*theta_ISA;
		
		return DISA;
	},
	
	OAT2DISA: function(OAT, delta){
		// Convert outside air temperature to delta ISA temperature
		// Units in Celsius
		//
		// Inputs:	OAT: outside air temperature (C)
		//			delta: pressure ratio
		// Outputs:	delta ISA temperature (C)

		var OAT0 = 288.15; // Kelvin
		var OAT = OAT + 273.15;
		var theta_ISA = this.Theta(delta, 0);
		var OAT_ISA = theta_ISA*OAT0;

		var DISA = OAT - OAT_ISA;
		
		return DISA;
	},

	Theta: function(delta, DISA){
		// Calculate temperature ratio (theta) given pressure ratio (delta) and 
		// delta ISA temperature in C
		//
		// Inputs:	delta: pressure ratio
		//			DISA: delta ISA temperature in C
		// Outputs:	temperature ratio (theta)

		var h_ft = this.PressureAltitude(delta, 0);
		var temp_ratio = 0;
		if (h_ft < 36089)
			temp_ratio = 1 - 6.87535e-6*h_ft + DISA/288.15;
		else
			temp_ratio = 0.7519 + DISA/288.15;
		
		return temp_ratio;
	},
	
	Delta: function(h_ft){
		// Calculate pressure ratio (delta) given pressure altitude in ft
		//
		// Inputs:	h_ft: pressure altitude in ft
		// Outputs:	pressure ratio (delta)

		var pressure_ratio = 0;
		if (h_ft < 36089)
			pressure_ratio = Math.pow((1 - h_ft/145442),5.2559);
		else
			pressure_ratio = 0.22336*Math.exp(-((h_ft-36089)/20806));
		
		return pressure_ratio;
	},
	
	Sigma: function(delta, theta){
		// Calculate the density ratio

		var sigma = delta/theta;
		
		return sigma;
	},
	
	Mu: function(theta, unit){
		// Calculate the dynamic viscosity of air
		//
		// Inputs:	theta: temperature ratio
		//			unit: 0 for imperial (slug/ft/sec), 1 for SI (Pa.sec)
		// Outputs:	dynamic viscosity (mu)

		var OAT0 = 288.15; // Kelvin
		var C1 = 1.458e-6; // kg/m/s/sqrt(K)
		var S = 110.4; // Kelvin
		var mu = Math.sqrt(OAT0)*C1*Math.pow(theta,(3/2))/(theta+S/OAT0);

		if (unit==0)
			mu = mu*0.0208854;
		
		return mu;
	},
	
	Nu: function(delta, theta, unit){
		// Calculate the kinematic viscosity of air
		//
		// Inputs:	delta: pressure ratio
		//			theta: temperature ratio
		//			unit: 0 for imperial (sq.ft/sec), 1 for SI (sq.m/sec)
		// Outputs:	kinematic viscosity (nu)

		var mu = this.Mu(theta, unit);
		var sigma = this.Sigma(delta, theta);
		var rho = this.Density(sigma, unit);
		var nu = mu/rho;
		
		return nu;
	},

	
	DynPres2TAS: function(q, sigma, unit){
		// Convert dynamic pressure to TAS
		//
		// Inputs:	q: dynamic pressure
		//			sigma: density ratio
		//			unit: 0 for imperial (kt, psf), 1 for SI (m/s, N/sq.m)
		// Outputs:	true airspeed (TAS)

		var rho = this.Density(sigma, unit);
		var TAS = Math.sqrt(2*q/rho);
		if (unit==0){
			// Convert to kts
			TAS = TAS*0.592484;
		}
		
		return TAS;
	},
	
	EAS2DynPres: function(EAS, unit){
		// Convert EAS to dynamic pressure
		//
		// Inputs:	EAS: equivalent airspeed
		//			unit: 0 for imperial (kt, psf), 1 for SI (m/s, N/sq.m)
		// Outputs:	dynamic pressure

		var rho0 = this.Density(1, unit);
		if (unit==0){
			// Convert to ft/s
			EAS = EAS/0.592484;
		}
		var q = 1/2*rho0*Math.pow(EAS,2);
		
		return q;
	},
	
	ImpactPres2CAS: function(qc, unit){
		// Convert impact pressure to CAS
		//
		// Inputs:	qc: impact pressure
		//			unit: 0 for imperial (kt, psf), 1 for SI (m/s, N/sq.m)
		// Outputs:	calibrated airspeed (CAS)

		var a0 = this.SpeedOfSound(1, unit);
		if (unit==0){
			// Convert to ft/s
			a0 = a0*1.68781;
		}
		var P0 = this.StaticPressure(1, unit);
		var CAS = a0*Math.sqrt(5*(Math.pow((qc/P0+1),(2/7))-1));
		if (unit==0){
			// Convert to kts
			CAS = CAS*0.592484;
		}
		
		return CAS;
	},
	
	Mach2DynPres: function(Mach, delta, unit){
		// Convert Mach to dynamic pressure
		//
		// Inputs:	Mach
		//			delta: pressure ratio
		//			unit: 0 for imperial (psf), 1 for SI (N/sq.m)
		// Outputs:	true airspeed (TAS)

		var P = this.StaticPressure(delta, unit);
		var q = 1/2*1.4*P*Math.pow(Mach,2);
		
		return q;
	},
	
	TAS2DynPres: function(TAS, sigma, unit){
		// Convert TAS to dynamic pressure
		//
		// Inputs:	TAS: true airspeed
		//			sigma: density ratio
		//			unit: 0=imperial (psf, kt), 1=SI (N/sq.m, m/s)
		// Outputs:	dynamic pressure

		var rho = this.Density(sigma, unit);
		if (unit==0){
			// Convert to ft/s
			TAS = TAS/0.592484;
		}
		var q = 1/2*rho*Math.pow(TAS,2);
		
		return q;
	},
	
	CAS2DynPres: function(CAS, delta, unit){
		// Convert CAS to dynamic pressure
		//
		// Inputs:	CAS: calibrated airspeed
		//			delta: pressure ratio
		//			unit: 0 for imperial (kt, psf), 1 for SI (m/s, N/sq.m)
		// Outputs:	dynamic pressure

		// Impact pressure
		var a0 = this.SpeedOfSound(1, unit);
		var P0 = this.StaticPressure(1, unit);
		if (unit==0){
			// Convert to ft/s
			CAS = CAS*1.68781;
			a0 = a0*1.68781;
		}
		var qc = (Math.pow((Math.pow((CAS/a0),2)/5+1),(7/2))-1)*P0;

		// P and Mach
		var P = this.StaticPressure(delta, unit);
		var Mach = Math.sqrt((Math.pow((qc/P+1),(2/7))-1)/0.2);

		q = 1/2*1.4*P*Math.pow(Mach,2);
		
		return q;
	},
	
	DynPres2EAS: function(q, unit){
		// Convert dynamic pressure to EAS
		//
		// Inputs:	q: dynamic pressure
		//			unit: 0 for imperial (kt, lbf), 1 for SI (m/s, N)
		// Outputs:	equivalent airspeed (EAS)

		var rho0 = this.Density(1, unit);
		var EAS = Math.sqrt(2*q/rho0);
		if (unit==0)
			EAS = EAS*0.592484;
		
		return EAS;
	},
	
	DynPres2ImpactPres: function(q, delta, unit){
		// Convert dynamic pressure to impact pressure
		//
		// Inputs:	q: dynamic pressure
		//			delta: pressure ratio
		//			unit: 0 for imperial (psf), 1 for SI (N/sq.m)
		// Outputs:	impact pressure

		var P = this.StaticPressure(delta, unit);
		var Mach = this.DynPres2Mach(q, delta, unit);
		var qc = P*(Math.pow((1+0.2*Math.pow(Mach,2)),(7/2))-1);
		
		return qc;
	},
	
	DynPres2Mach: function(q, delta, unit){
		// Convert dynamic pressure to Mach
		//
		// Inputs:	q: dynamic pressure
		//			delta: pressure ratio
		//			unit: 0 for imperial (psf), 1 for SI (N/sq.m)
		// Outputs:	Mach

		var P = this.StaticPressure(delta, unit);
		var Mach = Math.sqrt(2*q/1.4/P);
		
		return Mach;
	}
};