#include <TLorentzVector.h>

class kinematics {
	private:
		const double E1D_E0 = 4.802;
		const double MASS_E = 0.000511;
		const double MASS_P = 0.93827203;
		TVector3 e_mu_prime_3;
		TLorentzVector e_mu_prime;
		double _p;
		double _cx;
		double _cy;
		double _cz;

	public:
		kinematics() {
		}

		~kinematics() {
		}

		void SetVal(double p = 0 , double cx = 0 , double cy = 0 , double cz = 0){
			_p =  p; 
			_cx = cx; 
			_cy = cy; 
			_cz = cz; 

			e_mu_prime_3.SetXYZ(_p*_cx,p*_cy,p*_cz);	
			e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
		}

		double _Q2_calc(){
			TLorentzVector e_mu(0.0,0.0, sqrt(E1D_E0*E1D_E0-MASS_E*MASS_E), E1D_E0);
			TLorentzVector q_mu = (e_mu - e_mu_prime);
			return -q_mu.Mag2();
		}

		double _W_calc(){
			TLorentzVector e_mu(0.0,0.0, sqrt(E1D_E0*E1D_E0-MASS_E*MASS_E), E1D_E0);
			TLorentzVector q_mu = (e_mu - e_mu_prime);
			TVector3 p_mu_3(0,0,0);
			TLorentzVector p_mu;
			p_mu.SetVectM(p_mu_3,MASS_P);
			return (p_mu + q_mu).Mag();
		}

		double GetW() {
			return _W_calc();
		}
		double GetQ2() {
			return _Q2_calc();
		}

};