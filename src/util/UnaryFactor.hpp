#include<math.h>

//namespace graphchi{
const int ColorLevel = 12;
const int Len = ColorLevel*2+1;
//}
//	The table over a pair of variables
struct BinaryFactor {
	double *data;
	int n1, n2;

	BinaryFactor(int _n1 = 1, int _n2 = 1) : n1(_n1), n2(_n2) {
		data = new double[n1 * n2];
	}

	BinaryFactor(const BinaryFactor &other) {
		n1 = other.n1;
		n2 = other.n2;
		data = new double[n1 * n2];
		for (int i = 0; i < n1 * n2; i++)
			data[i] = other.data[i];
	}

	BinaryFactor &operator=(const BinaryFactor &other) {
		if (n1 * n2 != other.n1 * other.n2) {
			delete[] data;
			data = new double[other.n1 * other.n2];
		}
		n1 = other.n1;
		n2 = other.n2;

		for (int i = 0; i < n1 * n2; i++)
			data[i] = other.data[i];
		return *this;
	}

	void resize(int _n1, int _n2) {
		n1 = _n1;
		n2 = _n2;
		delete[] data;
		data = new double[n1 * n2];
	}

	inline void normalize() {
		double max_value = data[0];
		for (int i = 0; i < n1 * n2; i++)
			max_value = std::max(data[i], max_value);

		//	Scale and compute the normalized factor
		double Z = 0.0;
		for (int i = 0; i < n1 * n2; i++) {
			data[i] -= max_value;
			Z += exp(data[i]);
		}
		double logZ = log(Z);
		for (int i = 0; i < n1 * n2; i++)
			data[i] -= logZ;
	}
};

struct UnaryFactor{
	//double *data;	
	//double data[ColorLevel*2+1];	
	double data[Len];	
	//int n = Len;//ColorLevel*2+1;
	UnaryFactor(){
		/*
		n = 0;
		data = NULL;
		*/
		//n = Len;//ColorLevel*2+1;
	}
	UnaryFactor(int _n ){
		//n = Len;//ColorLevel*2+1;
		//data = new double[n];	
		//n = ColorLevel*2+1;	
	};
	
	void copy(UnaryFactor& other){
		//assert(n == other.n);
		for(int i=0; i<Len; i++){
			data[i] = other.data[i];	
		}
	}	
	/*	
	UnaryFactor& operator= (UnaryFactor& other){
		assert(n == other.n);
		for(int i=0; i<n; i++){
			data[i] = other.data[i];
		}		
		return *this;
	}		

	UnaryFactor& operator= (const UnaryFactor& other){
		std::cout<<"this.n"<<n<<"\t other.n"<<other.n<<std::endl;
		assert(n == other.n);
		for(int i=0; i<n; i++){
			data[i] = other.data[i];
		}		
		return *this;
	}		
	*/

	/*
	~UnaryFactor(){
		//delete[] data;			
	}*/		

	void resize(int _n){
		/*
		if(n != _n ){
			delete[] data;
			n = _n;
			data = new double[n];	
		}	
		*/
		
		if(_n != Len){	
			std::cout<<"this.n"<<Len<<" _n="<<_n<<std::endl;
			assert(false);
		}
	}		

	inline void uniform() {
		for (int i = 0; i < Len; i++)
			data[i] = 0.0;
	}

	inline void normalize() {
		double max_value = data[0];
		for (int i = 0; i < Len; i++) {
			max_value = std::max(max_value, data[i]);
		}

		double Z = 0.0;
		for (int i = 0; i < Len; i++) {
			data[i] -= max_value;
			Z += exp(data[i]);
		}
		double logZ = log(Z);
		for (int i = 0; i < Len; i++)
			data[i] -= logZ;
	}

	inline void times(const UnaryFactor& other) {
		for (int i = 0; i < Len; i++)
			data[i] += other.data[i];
	}

	inline void divide(const UnaryFactor& other) {
		for (int i = 0; i < Len; i++)
			data[i] -= other.data[i];
	}

	/** this(x) = sum_y fact(x,y) * other(y) */
	inline void convolve(const BinaryFactor& bin_fact, const UnaryFactor& other) {
		// Compute C(x) = Sum_y A(x,y) B(y)
		for (int i = 0; i < Len; i++) {
			double sum = 0.0f;
			//for (int j = 0; j < other.n; j++)
			for (int j = 0; j < Len; j++)
				sum += exp(bin_fact.data[i * bin_fact.n2 + j] + other.data[j]);

			data[i] = log(sum);
		}
	}

	inline void convolve(int uObs, int vObs, int deltaColor, int colors, double lambda, const UnaryFactor& other) {
		// Compute C(x) = Sum_y A(x,y) B(y)
		for (int i = 0; i < Len; i++) {
			double sum = 0.0f;
			//for (int j = 0; j < other.n; j++) {
			for (int j = 0; j < Len; j++) {
				double x1 = vObs + deltaColor * (i - colors);
				double x2 = uObs + deltaColor * (j - colors);
				double binFact = -fabs(x1 - x2) / lambda;

				sum += exp(binFact + other.data[j]);
			}
			data[i] = log(sum);
		}
	}

	inline double residual(const UnaryFactor& other) const {
		double residual = 0;
		for (int i = 0; i < Len; i++)
			//residual += fabs(data[i] - other.data[i]);
			residual = std::max(residual, fabs(data[i] - other.data[i]));

		return residual;
	}

	inline double norm() const {
		double norm = 0.0;
		double val = log(1.0 / Len);
		for (int i = 0; i < Len; i++)
			//norm += fabs(data[i] - val);
			norm = std::max(norm, fabs(data[i] - val));

		return norm;
	}

	/** This = other * damping + this * (1-damping) */
	inline void damp(const UnaryFactor& other, double damping) {
		if (damping == 0)
			return;
		for (int i = 0; i < Len; i++)
			data[i] = std::log(damping * std::exp(other.data[i]) + (1.0 - damping) * std::exp(data[i]));
	}

	/** get the max assignment */
	inline size_t max_asg() const {
		size_t max_asg = 0;
		double max_value = data[0];
		for (uint16_t i = 0; i < Len; i++) {
			if (data[i] > max_value) {
				max_value = data[i];
				max_asg = i;
			}
		}
		return max_asg;
	}		
};
