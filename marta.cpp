// Sexual selection and evolution of ornament
// COndition dependenrt project
// No kill anymore but reduced fecundity
// sempre by Marta


#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <random>
#include <chrono>
#include <iomanip>
#include <cassert>


std::chrono::high_resolution_clock::time_point tp =
std::chrono::high_resolution_clock::now();
unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());

// create and seed pseudo-random number generator
std::mt19937_64 rng(seed);


/**model parameters**/
int iNgenesQ = 10;
double m = 0.001; // mutational rate p and t
double ni = 0.005; // mutational rate di Q (deletirious)
double signma = 0.05; // distribution withd
double vi = 0.0001; // beneficial mutational rate di Q
double alfaA = 5.0; // effeiciency resource -> ornament
double alfaa = 4.0;
double viaA = 1.0;
double viaa = 0.5; // general viability 
double beta = 0.5;
double betaA = 0.5;
double yeta = 0.05;
int iPopsize = 200; // keep it even
int nTryMale = 100; // numeber of males female watch and compare
int tmax = 2001; // +1!!!
bool isfemale;
double f = 0.;
int iNsimulation = 1;

/*Definition of the class Female*/
class Female {
public:
	Female();
	std::vector <int> vectorQ;
	double getp() const { return (p); } // to caluclate average p 
	double getQ();
	std::pair <double, double> gett1t2() const { return std::pair <double, double>(t1, t2); } // to caluclate average p 
	void writegenome(double &, std::vector <double> &, std::vector <int> &);
private:
	double p; // first allele preference
	double t1; // 1 trait all
	double t2;
};
/*Definition of the class Male*/
class Male {
public:
	Male();
	std::vector <int> vectorQ;
	bool isAlive = true;
	double getp() const { return (p); } // to caluclate average p 
	double getQ();
	double getS();
	double gett();
	std::pair <double, double> gett1t2() const { return std::pair <double, double>(t1, t2); }
	void writegenome(double &, std::vector <double> &, std::vector <int> &);
private:
	double p;
	double t1;
	double t2;
};

/*Global variables*/



std::vector <Female> popFemale(iPopsize / 2);
std::vector <Male> popMale(iPopsize / 2);

std::vector <Female> popOffFem(iPopsize / 2);
std::vector <Male> popOffMale(iPopsize / 2);

std::vector <double> vecpref;
std::vector <double> vectrait1;
std::vector <double> vectrait2;
std::vector <double> vecTExpress;
std::vector <double> vecQualityfem;
std::vector <double> vecQualitymal;
std::vector <int> vecpopsizemale;
std::vector <int> vecpopsizefemale;


/*Implementation of the class Individual*/

Female::Female() {

	vectorQ = std::vector <int>(iNgenesQ, 0);
	std::normal_distribution <double> normalpre(0.0, 0.5);
	double y = normalpre(rng);
	double y2 = normalpre(rng);
	double y3 = normalpre(rng);
	p = y;
	t1 = y2;
	t2 = y3;

	// show(); // to show how preferences where set
}

Male::Male() {
	vectorQ = std::vector <int>(iNgenesQ, 0);
	std::normal_distribution <double> normalpre(0.0, 0.5);
	double y = normalpre(rng);
	double y2 = normalpre(rng);
	double y3 = normalpre(rng);
	p = y;
	t1 = y2;
	t2 = y3;
}

double Female::getQ() {
	double nQ = 0.;
	for (int x = 0; x < iNgenesQ; ++x) {
		double n = vectorQ[x];
		nQ += n;
	}
	double div = vectorQ.size();
	nQ /= div;
	return nQ;
}

double Male::getQ() {
	double nQ = 0.;
	for (int x = 0; x < iNgenesQ; ++x) {
		double n = vectorQ[x];
		nQ += n;
	}
	double div = vectorQ.size();
	nQ /= div;
	return nQ;
}

double Male::gett() {
	double qA = getQ();
	double t = (1 - f)*(qA*t1 + (1 - qA)*t2) + f*((t1 + t2) / 2);
	return t;
}

double Male::getS() {
	double qA = getQ();
	double t = gett();
	double dAlfa = alfaA * qA + alfaa* (1 - qA);
	double dS = t * dAlfa;
	return dS;
}

void Female::writegenome(double &P, std::vector <double> &vecT, std::vector <int>  &Nexvec)
{
	p = P;
	t1 = vecT[0];
	t2 = vecT[1];
	vectorQ = Nexvec;
	//std::cout << " Female DONE!and her p1 =" << p1 << "\n";
}


void Male::writegenome(double &P, std::vector <double> &vecT, std::vector <int>  &Nexvec)
{
	p = P;
	t1 = vecT[0];
	t2 = vecT[1];
	vectorQ = Nexvec;
	//std::cout << " Female DONE!and her p1 =" << p1 << "\n";
}

// Output files

std::ofstream oo("Individual.csv");
std::ofstream ofs("Alltrait.csv");


/*Simulation code*/
void variationquality(std::vector <Female> &pop) {
	for (int n = 0; n < (iNgenesQ*(iPopsize /3));) {
		std::uniform_int_distribution <int> ind(0, (iPopsize / 4) - 1);
		int i = ind(rng);
		std::uniform_int_distribution <int> row(0, iNgenesQ - 1);
		int x = row(rng);
		pop[i].vectorQ[x] = 1;
		++n;
	}

	for (int n = 0; n < (iNgenesQ*(iPopsize / 3)) ;) {
		std::uniform_int_distribution <int> ind((iPopsize / 4), (iPopsize / 2) - 1);
		int i = ind(rng);
		std::uniform_int_distribution <int> row(0, iNgenesQ - 1);
		int x = row(rng);
		pop[i].vectorQ[x] = 1;
		++n;
	}

}

void variationquality(std::vector <Male> &pop) {

	for (int n = 0; n < (iNgenesQ*(iPopsize / 3)) ;) {
		std::uniform_int_distribution <int> ind(0, (iPopsize / 4) - 1);

		int i = ind(rng);

		std::uniform_int_distribution <int> row(0, iNgenesQ - 1);

		int x = row(rng);
		pop[i].vectorQ[x] = 1;
		++n;
	}


	for (int n = 0; n < (iNgenesQ*(iPopsize / 3));) {
		std::uniform_int_distribution <int> ind((iPopsize / 4), (iPopsize / 2) - 1);
		int i = ind(rng);
		std::uniform_int_distribution <int> row(0, iNgenesQ - 1);
		int x = row(rng);

		pop[i].vectorQ[x] = 1;
		++n;
	}
}

void getaverageQ() {
	double sump = 0.;
	for (int i = 0; i < popFemale.size(); ++i) {
		double p = popFemale[i].getQ();
		sump += p;
	}
	double x = popFemale.size();
	sump /= x;
	//std::cout << " Average p in the female is"<< sump;
	vecQualityfem.push_back(sump); // keep track of it

	double sumq = 0.;
	for (int i = 0; i < popMale.size(); ++i) {
		double q = popMale[i].getQ();
		sumq += q;
	}
	double x1 = popMale.size();
	sumq /= x1;
	//std::cout << " Average p in the female is"<< sump;
	vecQualitymal.push_back(sumq); // keep track of i
}

void getaveragep() {
	double sump = 0.;
	for (int i = 0; i < popFemale.size(); ++i) {
		double p = popFemale[i].getp();
		sump += p;
	}
	double x = popFemale.size();
	sump /= x;
	//std::cout << " Average p in the female is"<< sump;
	vecpref.push_back(sump); // keep track of it
}

void getaveraget() {
	double sumt1 = 0.;
	double sumt2 = 0.;
	double sumtot = 0.;
	for (int i = 0; i < popMale.size(); ++i) {
		std::pair <double, double> t = popMale[i].gett1t2();
		double tot = popMale[i].gett();
		sumt1 += t.first;
		sumt2 += t.second;
		sumtot += tot;
	}
	double x = popMale.size();
	sumt1 /= x;
	sumt2 /= x;
	sumtot /= x;
	//std::cout << " Average t in the male is" << sumt << "\n";
	vectrait1.push_back(sumt1); // keep track of it
	vectrait2.push_back(sumt2);
	vecTExpress.push_back(sumtot);
}



void addmutation(double &val) { // mutation for alleles p and t
	std::bernoulli_distribution ismutation(m);
	if (ismutation(rng)) {
		std::normal_distribution <double> dis(0, signma);
		double mut = dis(rng);
		val += mut;
	}
}

void addQualMut(int &all) {
	std::bernoulli_distribution ismutation(ni);
	if (ismutation(rng)) {
		all = 0;
		//std::cout << "Bad quality MUTATION HAPPEN!";
	}

	std::bernoulli_distribution goodmut(vi);
	if (goodmut(rng)) {
		all = 1;
		//std::cout << "GOODDDDD quality MUTATION HAPPEN!";
	}
}

void madeoffspring(const Female &fem, const Male &male, const int &num, const bool &isfemale) {

		double newp; // write alleles to inherit
		std::vector <double> vecT(2, 0);
		std::vector <int>  vectorQ(iNgenesQ, 0);
		std::bernoulli_distribution Nallele(0.5); // mendelian inheritance

												  // select new alleles p
		double pfem = fem.getp();
		//std::cout << " \n p1 = " << pfem.first << " and p2 = " << pfem.second;
		double pmal = male.getp();

		int x = Nallele(rng);
		switch (x) {
		case 0:
			newp = pfem;
			break;
		case 1:
			newp = pmal;
			break;
		default:
			std::cout << " something wronk in chosing allele p for offspring";
		}

		addmutation(newp);

		//std::cout << " \n Alleles selected from mother" << vecP[0] << " And from the father = " << vecP[1] ;

		// select new alleles t
		std::pair <double, double> tfem = fem.gett1t2();
		std::pair <double, double> tmal = male.gett1t2();


		int y = Nallele(rng);
		switch (y) {
		case 0:
			vecT[0] = tfem.first;
			break;
		case 1:
			vecT[0] = tmal.first;
			break;
		default:
			std::cout << " something wronk in chosing allele t  for offspring";
		}

		addmutation(vecT[0]);

		int y1 = Nallele(rng);
		switch (y1) {
		case 0:
			vecT[1] = tfem.second;
			break;
		case 1:
			vecT[1] = tmal.second;
			break;
		default:
			std::cout << " something wronk in chosing allele t for offspring";
		}

		addmutation(vecT[1]);

		// select new quality matrix

		for (int i = 0; i < iNgenesQ; ++i) {
			int q = Nallele(rng);
			switch (q) {
			case 0:
				vectorQ[i] = fem.vectorQ[i];
				break;
			case 1:
				vectorQ[i] = male.vectorQ[i];
				break;
			default:
				std::cout << " something wronk in chosing allele t for offspring";
			}

			addQualMut(vectorQ[i]);

		}

		if (isfemale) {
			popOffFem[num].writegenome(newp, vecT, vectorQ); // first offspring is a female!! congrats
	
		}
		if (!isfemale) {
			popOffMale[num].writegenome(newp, vecT, vectorQ); // second offspring is a male!! congrats (but same father??, yes apparently she mate once)		
		
		}
}

void defineRates( std::vector <double> &rates ) {
	
	for (int i = 0; i < popMale.size(); ++i) {
		double qA = popMale[i].getQ();
		double t = popMale[i].gett();
		double v = viaA*qA + viaa* (1 - qA);
		double b = betaA*qA + beta* (1 - qA);
		double cm = 1 - exp(-b*(t*t));
		double hm = v * (1 - cm);
		rates.push_back(hm);
	}

}

int selectmale(const Female &fem, const std::vector <double> &rates) {
	//std::cout << " \n for this female is :";
	
	std::vector <int> nMale;
	std::vector <double> attract;
	// selec nTryMale and write their attractivness


	for (int i = 0; i < nTryMale; ++i) {
		std::discrete_distribution <int> Lottery(rates.begin(), rates.end());
		int next = Lottery(rng);
		nMale.push_back(next);
		// calculate attractivness
		double p = fem.getp();
		double s = popMale[next].getS();
		double r = exp(p*s);
		attract.push_back(r);

	}

	std::discrete_distribution <int> weightLottery(attract.begin(), attract.end());
	int j = weightLottery(rng);

	return nMale[j];

}

int selectfemale() {

	std::vector <double> rates;
	// calculate fecundity
	for (int i = 0; i < popFemale.size(); ++i) {
		double qA = popFemale[i].getQ();
		double p = popFemale[i].getp();
		double v = viaA*qA + viaa* (1 - qA);
		double cf = 1 - exp(-yeta*(p*p));
		double hf = v * (1 - cf);
		rates.push_back(hf);
	}

	std::discrete_distribution <int> weightLottery(rates.begin(), rates.end());
	int j = weightLottery(rng);
	return j;

}

void NextGeneration() {
	std::vector <double> rates;
	defineRates(rates);


	isfemale = true;
	for (int i = 0; i < popOffFem.size(); ++i) {
		int f = selectfemale(); // selected based on fecundity
		int y = selectmale(popFemale[f], rates); // number of males will mate with female f
		madeoffspring(popFemale[f], popMale[y], i, isfemale);
	}
	isfemale = false;

	for (int i = 0; i < popOffMale.size(); ++i) {
		int f1 = selectfemale(); // selected based on fecundity
		int y1 = selectmale(popFemale[f1], rates); // number of males will mate with female f
		madeoffspring(popFemale[f1], popMale[y1], i, isfemale);
	}

}



void update() {
	popFemale = popOffFem;
	popMale = popOffMale;
	popOffFem = std::vector <Female>(iPopsize/2);
	popOffMale = std::vector <Male>(iPopsize/2);

	getaveragep();
	getaveraget();
	getaverageQ();
	vecpopsizemale.push_back(popMale.size());
	vecpopsizefemale.push_back(popFemale.size());
}



void writeshit() {

	oo << "Male \n";
	oo << "Q" << ',' << "T1" << ',' << "T2" << ',' << "T" << ',' << std::endl;
	for (int i = 0; i < popMale.size(); ++i) {
		std::pair <double, double>  t = popMale[i].gett1t2();
		double s = popMale[i].gett();
		double Q = popMale[i].getQ();
		oo << Q << ',' << t.first << ',' << t.second << ',' << s << ',' << std::endl;
	}
	oo << "Female  \n";
	oo << "q" << ',' << "p" << std::endl;
	for (int i = 0; i < popFemale.size(); ++i) {
		double p = popFemale[i].getp();
		double Q = popFemale[i].getQ();
		oo << Q << ',' << p << ',' << std::endl;
	}

}

void refresh() {


	popFemale = std::vector <Female>(iPopsize / 2);
	popMale = std::vector <Male>(iPopsize / 2);

	assert(popFemale.size() == iPopsize / 2);
	assert(popMale.size() == iPopsize / 2);

	popOffFem = std::vector <Female>(iPopsize / 2);
	popOffMale = std::vector <Male>(iPopsize / 2);

	assert(popOffFem.size() == iPopsize / 2);
	assert(popOffMale.size() == iPopsize / 2);


	vecpref.erase(vecpref.begin(), vecpref.end());
	assert(vecpref.empty());

	vectrait1.erase(vectrait1.begin(), vectrait1.end());
	assert(vectrait1.empty());

	vectrait2.erase(vectrait2.begin(), vectrait2.end());
	assert(vectrait2.empty());

	vecTExpress.erase(vecTExpress.begin(), vecTExpress.end());
	assert(vecTExpress.empty());

	vecQualityfem.erase(vecQualityfem.begin(), vecQualityfem.end());
	assert(vecQualityfem.empty());

	vecQualitymal.erase(vecQualitymal.begin(), vecQualitymal.end());
	assert(vecQualitymal.empty());

	vecpopsizemale.erase(vecpopsizemale.begin(), vecpopsizemale.end());
	assert(vecpopsizemale.empty());

	vecpopsizefemale.erase(vecpopsizefemale.begin(), vecpopsizefemale.end());
	assert(vecpopsizefemale.empty());

}
/*Main (se mai ci arriverò)(sono già qua) ( e siamo a 350 lines!)(450 now baby)*/

int main() {


	// obtain seed from system clock

	for (int ns = 0; ns < iNsimulation; ++ns) {
		// initial shit
		std::cout << " SIMULATION N =  " << ns << std::endl;
		oo << " \n Simulation N = " << ns << "\n";
		if (ns != 0)
			refresh();

		variationquality(popMale); // distributred high quality genes around males
								   //std::cout << " FIRST VARIATION QUALITY WORKED \n";
		variationquality(popFemale); //same for females
									 //showmatrixQ(); // check it works, it does

	
									 //calcualte average p in females
		getaveragep();
		getaveraget();
		getaverageQ();
		vecpopsizemale.push_back(popMale.size());
		vecpopsizefemale.push_back(popFemale.size());
	

		/*START SIMULATION FOR A LOT OF GENERATIONSSSSSS*/

		for (int t = 0; t < tmax; ++t) {


			if (t % 100 == 0) {
				std::cout << " generation N =  " << t << std::endl;
				oo << " \n NEXT GENERATION " << t << "\n";
				writeshit();
				//std::cout << " POp female = " << popFemale.size() << " POp male = " << popMale.size() << std::endl;
			}

			if (popFemale.size() == 0 || popMale.size() < nTryMale) {
				std::cout << " POPULATION GOT EXTINTTT MAAAAAAAANNN" << std::endl;
				break;
			}

			// creat next generation
			NextGeneration();
	
			// Survivors became PARENTS!! congratss and CLEAR OFF POP, save p and t traits
			update();


		}
	

		// Check evolution of p and t
		ofs << " SIMULATION N = " << ns << " Random seed = " << seed << std::endl;

		ofs << " Time t " << ',' << " Preference p " << ',' << " Trait t1  " << ',' << " Trait t2 " << ',' << " t express " << ',' << " Population size male  " << ',' << " Population size FEMmale  " << ',' << " Mean quality female " << ',' << " Mean quality male " << std::endl;
		for (int i = 0; i < tmax + 1; ) {
			//std::cout << i << ',' << vecpref[i] << ',' << vectrait1[i] << ',' << vectrait2[i] << ',' << vecpopsizemale[i] << ',' << vecpopsizefemale[i] << ',' << vecQualityfem[i] << ',' << vecQualitymal[i] << std::endl;
			ofs << i << ',' << vecpref[i] << ',' << vectrait1[i] << ',' << vectrait2[i] << ',' << vecTExpress[i] << ',' << vecpopsizemale[i] << ',' << vecpopsizefemale[i] << ',' << vecQualityfem[i] << ',' << vecQualitymal[i] << std::endl;
			i = i + 100;
		}

}
	return 0;
}

