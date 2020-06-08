package model;

import com.jsyn.ports.UnitOutputPort;
import com.jsyn.unitgen.Circuit;

public class Clavinet extends Circuit {
	
	private long n = 0;
	
	public UnitOutputPort output1;
	public UnitOutputPort output2;
	
	private double fs = 44100.0;
	private double k = 1 / fs;
	private double dur; // note duration
	private double coeffDur = 1;
	
	private double rho; // density [kg/m^3]
	private double T; // tension [N]
	private double radius; // cross-section radius [m]
	private double A; // cross-section area [m^2]
	private double I; // moment of inertia
	private double L; // string length [m]
	private double E; // Young's modulus
	private double K; // stiffness parameter
	private double c; // wave speed
	private double sigma1; // frequency dependent damping
	private int N; // number of space steps
	private double h; // space step size
	
	private double yarnPos; // percentage of yarn position
	
	private final int maxLength = 10000; // max number of space steps
	private double sY; // yarn damping
	private double sS; // string damping
	private double[] sigma0; // frequency independent damping
	
	private double[] uNext;
	private double[] u;
	private double[] uPrev;
	
	private double exciteLength; // portion of string excited
	
	private double[] J;
	
	private int outPos1 = 0;
	private int outPos2 = 0;
	
	private double[] coeff0;
	private double[] coeff1;
	private double[] coeff2;
	private double[] coeff3;
	private double[] coeff4;
	
	private double tSilence1;
	private double tAttack;
	private double tRelease;
	private double tSilence2;
	private double force;
	
	private double epsilon;
	private double kStud;
	private double alpha;
	
	private double limiterValue;
	private boolean distortionActive;
	

	public Clavinet() {
		setupOutputPorts();
		setDuration(5.0);
		setDensity(8050.0);
		setTension(80.0);
		setRadius(0.0003556);
		setYoungModulus(174000000.0);
		setArea();
		setInertiaMoment();
		setWaveSpeed();
		setStiffnessParameter();
		setFrequencyDependentDamping(0.012);
		setYarnPercentagePosition(10.0);
		setYarnDamping(515.0);
		setStringDamping(0.005);
		stabilityCondition();
		setStringLength(1.0);
		updateSigma0();
		
		setExciteLength(0.7);
		settSilence1(0.0);
		settAttack(0.001 * fs);
		settRelease(0.001 * fs);
		settSilence2(0.0);
		setForce(50.0);
		setEpsilon(0.001);
		setkStud(10000);
		setAlpha(1.1);
		
		setLimiterValue(0.001);
		setDistortionActive(false);
		
		
		initializeUNext();
		initializeU();
		initializeUPrev();
		updateJ();
		updateOutPos();
		updateCoeff0();
		updateCoeff1();
		updateCoeff2();
		updateCoeff3();
		updateCoeff4();
		
		//System.out.println(rho+" "+T+" "+radius+" "+A+" "+I+" "+L+" "+E+" "+K+" "+c+" "+sigma1+" "+N+" "+h+" "+yarnPos+" "+sY+" "+sS);
		/*double sum = 0;
		for(int i=0; i<N; i++) sum += sigma0[i];
		System.out.println(sum);*/
	}
	
	public synchronized void setDuration(double durationSeconds) {
		dur = durationSeconds * fs;
	}
	
	public double getDuration() {
		return dur / fs;
	}
	
	public double getNumberOfSpaceSteps() {
		return N;
	}
	
	public double getSpaceStep() {
		return h;
	}
	
	public double getSampleFrequency() {
		return fs;
	}
	
	public double getTimeStep() {
		return k;
	}
	
	public double[] getUNext() {
		return uNext;
	}
	
	public double[] getU() {
		return u;
	}
	
	public double[] uPrev() {
		return uPrev;
	}
	
	public double getDensity() {
		return rho;
	}
	
	public synchronized void setDensity(double rho) {
		this.rho = rho;
		setStiffnessParameter();
		setWaveSpeed();
		stabilityCondition();
		updateJ();
	}
	
	public double getTension() {
		return T;
	}
	
	public synchronized void setTension(double tension) {
		this.T = tension;
		setWaveSpeed();
	}
	
	public double getRadius() {
		return radius;
	}
	
	public synchronized void setRadius(double radius) {
		this.radius = radius;
		setArea();
		setInertiaMoment();
	}
	
	public double getArea() {
		return A;
	}
	
	private synchronized void setArea() {
		A = Math.PI * radius * radius;
		setStiffnessParameter();
		setWaveSpeed();
		stabilityCondition();
		updateJ();
	}
	
	public double getInertiaMoment() {
		return I;
	}
	
	private synchronized void setInertiaMoment() {
		I = (Math.PI * Math.pow(radius, 4.0)) / 4.0;
		setStiffnessParameter();
		stabilityCondition();
		
	}
	
	public double getStringLength() {
		return L;
	}
	
	public synchronized void setStringLength(double stringLength) {
		this.L = stringLength;
		stabilityCondition();
		updateJ();
	}
	
	public double getYoungModulus() {
		return E;
	}
	
	public synchronized void setYoungModulus(double youngModulus) {
		this.E = youngModulus;
		setStiffnessParameter();
		stabilityCondition();
	}
	
	public double getStiffnessParameter() {
		return K;
	}
	
	private synchronized void setStiffnessParameter() {
		K = Math.sqrt(E*I/rho/A);
		updateCoeff4();
	}
	
	public double getWaveSpeed() {
		return c;
	}
	
	private synchronized void setWaveSpeed() {
		c = Math.sqrt(T/rho/A);
		stabilityCondition();
		updateCoeff2();
	}
	
	public double getFrequencyDependentDamping() {
		return sigma1;
	}
	
	public synchronized void setFrequencyDependentDamping(double sigma1) {
		this.sigma1 = sigma1;
		stabilityCondition();
		updateCoeff2();
		updateCoeff3();
	}
	
	public double getYarnPercentagePosition() {
		return yarnPos;
	}
	
	public synchronized void setYarnPercentagePosition(double yarnPos) {
		this.yarnPos = yarnPos;
		updateSigma0();
		updateOutPos();
	}
	
	public double getYarnDamping() {
		return sY;
	}
	
	public synchronized void setYarnDamping(double sY) {
		this.sY = sY;
		updateSigma0();
	}
	
	public double getStringDamping() {
		return sS;
	}
	
	public synchronized void setStringDamping(double sS) {
		this.sS = sS;
		updateSigma0();
	}
	
	public double[] getFrequencyIndependentDamping() {
		return sigma0;
	}
	
	private synchronized void updateSigma0() {
		sigma0 = new double[maxLength];
		int yarn = (int)(N * (yarnPos/100.0));
		for(int i=0; i<N; i++) {
			if(i < yarn) sigma0[i] = sY;
			else sigma0[i] = sS;
		}
		
		updateJ();
		updateCoeff0();
		updateCoeff1();
		updateCoeff2();
		updateCoeff3();
		updateCoeff4();
	}
	
	public double getExciteLength() {
		return exciteLength;
	}
	
	public synchronized void setExciteLength(double percentage) {
		exciteLength = L * percentage;
		updateJ();
	}
	
	public double[] getJ() {
		return J;
	}
	
	private synchronized void updateJ() {
		double excitePos = (L - exciteLength) / h;
		int excitePosFloor = (int) excitePos;
		double beta = excitePos - excitePosFloor;
		double[] delta = new double[maxLength];
		for(int i=0; i<N; i++) {
			if(i == excitePosFloor) delta[i] = (1-beta) / h;
			else if(i == excitePosFloor+1) delta[i] = beta / h;
			else delta[i] = 0;
		}
		J = new double[maxLength];
		for(int i=0; i<N; i++) {
			J[i] = k*k*delta[i]/(1-sigma0[i]*k)/(rho*A);
		}
	}
	
	public double[] getCoeff0() {
		return coeff0;
	}
	
	private synchronized void updateCoeff0() {
		coeff0 = new double[maxLength];
		for(int i=0; i<N; i++) {
			coeff0[i] = 2/(1+k*sigma0[i]);
		}
	}
	
	public double[] getCoeff1() {
		return coeff1;
	}
	
	private synchronized void updateCoeff1() {
		coeff1 = new double[maxLength];
		for(int i=0; i<N; i++) {
			coeff1[i] = (k*sigma0[i] - 1)/(k*sigma0[i] + 1);
		}
	}
	
	public double[] getCoeff2() {
		return coeff2;
	}
	
	private synchronized void updateCoeff2() {
		coeff2 = new double[maxLength];
		for(int i=0; i<N; i++) {
			coeff2[i] = (c*c*k*k + 2*sigma1*k)/((h*h)*(1+k*sigma0[i]));
		}
	}
	
	public double[] getCoeff3() {
		return coeff3;
	}
	
	private synchronized void updateCoeff3() {
		coeff3 = new double[maxLength];
		for(int i=0; i<N; i++) {
			coeff3[i] = 2*sigma1*k/((h*h)*(1+k*sigma0[i]));
		}
	}
	
	public double[] getCoeff4() {
		return coeff1;
	}
	
	private synchronized void updateCoeff4() {
		coeff4 = new double[maxLength];
		for(int i=0; i<N; i++) {
			coeff4[i] = (K*K*k*k)/((h*h*h*h)*(1+k*sigma0[i]));
		}
	}
	
	public double gettSilence1() {
		return tSilence1;
	}

	public synchronized void settSilence1(double tSilence1) {
		this.tSilence1 = tSilence1;
	}

	public double gettAttack() {
		return tAttack;
	}

	public synchronized void settAttack(double tAttack) {
		this.tAttack = tAttack;
	}

	public double gettRelease() {
		return tRelease;
	}

	public synchronized void settRelease(double tRelease) {
		this.tRelease = tRelease;
	}

	public double gettSilence2() {
		return tSilence2;
	}

	public synchronized void settSilence2(double tSilence2) {
		this.tSilence2 = tSilence2;
	}

	public double getForce() {
		return force;
	}

	public synchronized void setForce(double force) {
		this.force = force;
	}

	public double getEpsilon() {
		return epsilon;
	}

	public synchronized void setEpsilon(double epsilon) {
		this.epsilon = epsilon;
	}

	public double getkStud() {
		return kStud;
	}

	public synchronized void setkStud(double kStud) {
		this.kStud = kStud;
	}

	public double getAlpha() {
		return alpha;
	}

	public synchronized void setAlpha(double alpha) {
		this.alpha = alpha;
	}
	
	public double getLimiterValue() {
		return this.limiterValue;
	}
	
	public synchronized void setLimiterValue(double limiterValue) {
		this.limiterValue = limiterValue;
	}
	
	public boolean isDistortionActive() {
		return this.distortionActive;
	}
	
	public void setDistortionActive(boolean distortionActive) {
		this.distortionActive = distortionActive;
	}
	
	private synchronized void setupOutputPorts() {
		output1 = new UnitOutputPort();
		this.addPort(output1);
		
		output2 = new UnitOutputPort();
		this.addPort(output2);
	}
	
	private synchronized void stabilityCondition() {
		double coeff = c*c*k*k + 4*sigma1*k;
		h = Math.sqrt((coeff + Math.sqrt(coeff*coeff+((16*E*I*k*k)/(rho*A))))/2.0)*1.2;
		N = (int)(L/h);
	}
	
	
	
	private synchronized void initializeUNext() {
		uNext = new double[maxLength];
		for(int i=0; i<N; i++) uNext[i] = 0;
		stabilityCondition();
	}
	
	private synchronized void initializeU() {
		u = new double[maxLength];
		for(int i=0; i<N; i++) u[i] = 0;
		stabilityCondition();
	}
	
	private synchronized void initializeUPrev() {
		uPrev = new double[maxLength];
		for(int i=0; i<N; i++) uPrev[i] = 0;
		stabilityCondition();
	}
	
	private double forceInput(long t) {
		double tSustain = dur - tSilence1 - tSilence2 - tAttack - tRelease;
		double t1 = tSilence1+tAttack;
		double t2 = tSilence1+tAttack+tSustain;
		double t3 = tSilence1+tAttack+tSustain+tRelease;
		
		if(t > tSilence1 && t < t1) return 0.5*(1-Math.cos(Math.PI*(t-tSilence1)/tAttack))*force;
		else if(t >= t1 && t <= t2) return force;
		else if(t >= t2 && t2 <= t3) return 0.5*(1-Math.cos(Math.PI*(t-t3)/tRelease))*force;
		else return 0;
	}
	
	
	
	private double forceStud() {
		double excitePos = (L - exciteLength) / h;
		int excitePosFloor = (int) excitePos;
		double beta = excitePos - excitePosFloor;
		double studParam = -(1-beta)*u[excitePosFloor] - beta*u[excitePosFloor+1] - epsilon;
		if(studParam > 0) return kStud*Math.pow(studParam, alpha);
		else return 0;
	}
	
	private synchronized void updateOutPos() {
		outPos1 = (int)((N - N*(yarnPos/100.0)) / 2.0);
		outPos2 = (int)( (N*(yarnPos/100.0)) / 2.0 );
	}
	
	
	
	private synchronized void updateUNext() {
		double fIn = forceInput(n);
		double fStud = forceStud();
		for(int l=2; l<N-2; l++) {
			double d4 = u[l+2] - 4*u[l+1] + 6*u[l] - 4*u[l-1] + u[l-2];
			uNext[l] = coeff0[l]*u[l] + coeff1[l]*uPrev[l] + coeff2[l]*(u[l+1] - 2 * u[l] + u[l-1]) - coeff3[l]*(uPrev[l+1] - 2*uPrev[l] + uPrev[l-1]) - coeff4[l]*d4 + J[l]*(fStud - fIn);
		}
	}
	
	private synchronized void updateU() {
		for(int i=0; i<N; i++) {
			u[i] = uNext[i];
		}
	}
	
	private synchronized void updateUPrev() {
		for(int i=0; i<N; i++) {
			uPrev[i] = u[i];
		}
	}
	
	private synchronized void updateCoeffDur() {
		if(n >= dur) coeffDur = 0;
		else {
			if(dur-n > 10000) coeffDur = 50;
			else coeffDur = 50*(dur-n)/10000;
		}
	}
	
	public synchronized void excite() {
		/*double sum = 0;
		for(int i=0; i<N; i++) sum += uNext[i];
		System.out.println(sum);*/
		
		initializeUNext();
		initializeU();
		initializeUPrev();
		updateJ();
		n = 0;
		//updateCoeffDur();
	}
	
	public synchronized void deexcite() {
		n = (long)(2*dur);
	}
	

	@Override
	public synchronized void generate(int start, int end) {
		double[] out1 = output1.getValues();
		double[] out2 = output2.getValues();
		for(int i=start; i<end; i++) {
			updateUNext();
			updateUPrev();
			updateU();
			updateCoeffDur();
			//System.out.println(uNext[i]);
			out1[i] = coeffDur * uNext[outPos1];
			out2[i] = coeffDur * uNext[outPos2];
			if(this.distortionActive) {
				out1 = distortion(out1);
				out2 = distortion(out2);
			}
			n++;
		}
		
	}
	
	private double[] distortion(double[] vector) {
		for(int i=0; i<vector.length; i++) {
			if(vector[i] > this.limiterValue) {
				vector[i] = this.limiterValue;
			}
			else if(vector[i] < -this.limiterValue) {
				vector[i] = -this.limiterValue;
			}
		}
		return vector;
	}
	
	
	
}
