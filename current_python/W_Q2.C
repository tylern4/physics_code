#include <TLorentzVector.h>

class W_Q2 {
	private:
		double E1D_E0 = 4.802;
		double MASS_E = 0.000511;
		double MASS_P = 0.93827203;
		TVector3 e_mu_prime_3;
		TLorentzVector e_mu_prime;
		double _p;
		double _cx;
		double _cy;
		double _cz;
		double _W;
		double _Q2;

	public:
		W_Q2 () {
		}

		void Q2_calc(){
			TLorentzVector e_mu(0.0,0.0, sqrt(E1D_E0*E1D_E0-MASS_E*MASS_E), E1D_E0);
			TLorentzVector q_mu = (e_mu - e_mu_prime);
			_Q2 = -q_mu.Mag2();
		}

		void W_calc(){
			TLorentzVector e_mu(0.0,0.0, sqrt(E1D_E0*E1D_E0-MASS_E*MASS_E), E1D_E0);
			TLorentzVector q_mu = (e_mu - e_mu_prime);
			TVector3 p_mu_3(0,0,0);
			TLorentzVector p_mu;
			p_mu.SetVectM(p_mu_3,MASS_P);
			_W = (p_mu + q_mu).Mag();
		}

		void SetVal(double p = 0 , double cx = 0 , double cy = 0 , double cz = 0){
			_p =  p; 
			_cx = cx; 
			_cy = cy; 
			_cz = cz; 

			e_mu_prime_3.SetXYZ(_p*_cx,p*_cy,p*_cz);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

			W_calc();
			Q2_calc();

		}

		double GetW() {
			return _W;
		}
		double GetQ2() {
			return _Q2;
		}


};