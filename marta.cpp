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
double m = 0.01; // mutational rate p and t
double ni = 0.05; // mutational rate di Q (deletirious)
double signma = 0.05; // distribution withd
double alfaA = 5.0; // effeiciency resource -> ornament
double alfaa = 4.0;
double viaA = 1.0;
double viaa = 0.5; // general viability 
double beta = 0.5;
double betaA = 0.5;
double yeta = 0.05;
int iPopsize = 200; // keep it even
int tmax = 2001; // +1!!!
bool isfemale;
double f = 0.; // f for condition dependent expression of t 1= no condition dependent
double e = 0; // for environmental dependency on quality 0= no environmental effects
double emin = 0.0;
double emax = 1.0;
double vartp = 0.5; // witdh initial variation for genes trait p and t
int iNsimulation = 1;

/*Definition of the class Female*/
class Female {
public:
	Female();
	double getp() const { return (p); } // to caluclate average p 
	double getnumberQ() const { return (nQ); }
	double getQE() const { return (nQE); }
	double getQ();
	std::pair <double, double> gett1t2() const { return std::pair <double, double>(t1, t2); } // to caluclate average p 
	void writegenome(double &, std::vector <double> &, int  &);
private:
	double p; // first allele preference
	double t1; // 1 trait all
	double t2;
	double nQ;
	double nQE; // environmental quality
};
/*Definition of the class Male*/
class Male {
public:
	Male();
	double getp() const { return (p); } // to caluclate average p 
	double getnumberQ() const { return (nQ); }
	double getQE() const { return (nQE); }
	double getQ();
	double getS();
	double gett();
	std::pair <double, double> gett1t2() const { return std::pair <double, double>(t1, t2); }
	void writegenome(double &, std::vector <double> &, int  &);
private:
	double p;
	double t1;
	double t2;
	double nQ;
	double nQE; // environmental quality
};

/*Global variables*/


std::vector <Female> popFemale(iPopsize / 2);
std::vector <Male> popMale(iPopsize / 2);
std::vector <Female> popOffFem(iPopsize / 2);
std::vector <Male> popOffMale(iPopsize / 2);



/*Implementation of the class Individual*/

Female::Female() {

	std::uniform_int_distribution <int> qdist(iNgenesQ / 2, iNgenesQ);
	std::normal_distribution <double> normalpre(0.0, vartp);
	std::uniform_real_distribution <double> enviDist(emin, emax);
	p = normalpre(rng);
	t1 = normalpre(rng);
	t2 = normalpre(rng);
	nQ = qdist(rng);
	nQE = enviDist(rng);
}

Male::Male() {
	std::uniform_int_distribution <int> qdist(iNgenesQ / 2, iNgenesQ);
	std::normal_distribution <double> normalpre(0.0, vartp);
	std::uniform_real_distribution <double> enviDist(emin, emax);
	p = normalpre(rng);
	t1 = normalpre(rng);
	t2 = normalpre(rng);
	nQ = qdist(rng);
	nQE = enviDist(rng);
}

double Female::getQ() {
	double nQ = getnumberQ();
	double fq = nQ / iNgenesQ; // fraction good quality genes

	double QE = getQE();
	double quality = (1 - e)*(fq)+e*QE;

	return quality;
}

double Male::getQ() {
	double nQ = getnumberQ();
	double fq = nQ / iNgenesQ; 
	double QE = getQE();
	double quality = (1 - e)*(fq)+e*QE;
	return quality;
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

void Female::writegenome(double &P, std::vector <double> &vecT, int  &Nexqua)
{
	std::uniform_real_distribution <double> enviDist(emin, emax);
	p = P;
	t1 = vecT[0];
	t2 = vecT[1];
	nQ = Nexqua;
	nQE = enviDist(rng);
}

void Male::writegenome(double &P, std::vector <double> &vecT, int  &Nexqua)
{
	std::uniform_real_distribution <double> enviDist(emin, emax);
	p = P;
	t1 = vecT[0];
	t2 = vecT[1];
	nQ = Nexqua;
	nQE = enviDist(rng);
}


// Output files

std::ofstream oo;
std::ofstream ofs;
std::ofstream para;

/*Simulation code*/

double getaverageQ() {
	double sump = 0.;
	for (int i = 0; i < popFemale.size(); ++i) {
		double p = popFemale[i].getQ();
		sump += p;
	}
	double x = popFemale.size();
	sump /= x;

	double sumq = 0.;
	for (int i = 0; i < popMale.size(); ++i) {
		double q = popMale[i].getQ();
		sumq += q;
	}

	double x1 = popMale.size();
	sumq /= x1;
	double  tot = (sump + sumq) / 2;
	return tot;
}

double getaveragep() {
	double sump = 0.;
	for (int i = 0; i < popFemale.size(); ++i) {
		double p = popFemale[i].getp();
		sump += p;
	}
	double x = popFemale.size();
	sump /= x;
	return sump; // keep track of it
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
	ofs << ',' << sumt1 << ',' << sumt2 << ',' << sumtot;
}



void addmutation(double &val) { // mutation for alleles p and t
	std::bernoulli_distribution ismutation(m);
	if (ismutation(rng)) {
		std::normal_distribution <double> dis(0, signma);
		double mut = dis(rng);
		val += mut;
	}
}



void madeoffspring(Female &fem, Male &male, const int &num, const bool &isfemale) {

	double newp; // write alleles to inherit
	std::vector <double> vecT(2, 0);
	std::bernoulli_distribution Nallele(0.5); // mendelian inheritance

											  // select new alleles p
	double pfem = fem.getp();
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

	double pm = fem.getQ();
	double pf = male.getQ();


	std::binomial_distribution <int> nMotherQuality(iNgenesQ / 2, pm);
	std::binomial_distribution <int> nFatherQuality(iNgenesQ / 2, pf);

	int qm = nMotherQuality(rng);

	int qf = nFatherQuality(rng);


	std::poisson_distribution <int> badmutation(ni);
	double mut = badmutation(rng);


	int newq = qm + qf - mut;


	if (isfemale) {
		popOffFem[num].writegenome(newp, vecT, newq); // first offspring is a female!! congrats

	}
	if (!isfemale) {
		popOffMale[num].writegenome(newp, vecT, newq); // second offspring is a male!! congrats (but same father??, yes apparently she mate once)		

	}
}

void defineRates(std::vector <double> &rates) {

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

void defineFemaleFecundity(std::vector <double> &rates) {

	for (int i = 0; i < popFemale.size(); ++i) {
		double qA = popFemale[i].getQ();
		double p = popFemale[i].getp();
		double v = viaA*qA + viaa* (1 - qA);
		double cf = 1 - exp(-yeta*(p*p));
		double hf = v * (1 - cf);
		rates[i] = hf;
	}

}


int selectmale(const Female &fem, const std::vector <double> &rates) {

	std::vector <double> attract(popMale.size());
	// selec nTryMale and write their attractivness


	for (int i = 0; i < popMale.size(); ++i) {
		// calculate attractivness
		double p = fem.getp();
		double s = popMale[i].getS();
		double r = exp(p*s);
		attract[i] = r * rates[i];
	}

	std::discrete_distribution <int> Lottery(attract.begin(), attract.end());
	int j = Lottery(rng);
	return j;

}

int selectfemale(const std::vector <double> &rates) {

	std::discrete_distribution <int> weightLottery(rates.begin(), rates.end());
	int j = weightLottery(rng);
	return j;

}
void defineMaleFecundity(std::vector <double> &rates) {

	for (int i = 0; i < popMale.size(); ++i) {
		double qA = popMale[i].getQ();
		double t = popMale[i].gett();
		double v = viaA*qA + viaa* (1 - qA);
		double b = betaA*qA + beta* (1 - qA);
		double cm = 1 - exp(-b*(t*t));
		double hm = v * (1 - cm);
		rates[i] = hm;
	}

}

void NextGeneration() {
	std::vector <double> MaleRates(popMale.size());
	defineMaleFecundity(MaleRates);

	std::vector <double> FemaleRates(popFemale.size());
	defineFemaleFecundity(FemaleRates);

	isfemale = true;

	for (int i = 0; i < popOffFem.size(); ++i) {
		int f = selectfemale(FemaleRates); // selected based on fecundity
		int y = selectmale(popFemale[f], MaleRates); // number of males will mate with female f
		madeoffspring(popFemale[f], popMale[y], i, isfemale);
	}


	isfemale = false;

	for (int i = 0; i < popOffFem.size(); ++i) {
		int f = selectfemale(FemaleRates); // selected based on fecundity
		int y = selectmale(popFemale[f], MaleRates); // number of males will mate with female f
		madeoffspring(popFemale[f], popMale[y], i, isfemale);
	}
}



void update() {
	popFemale = popOffFem;
	popMale = popOffMale;
	popOffFem = std::vector <Female>(iPopsize / 2);
	popOffMale = std::vector <Male>(iPopsize / 2);
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

void writeState() {
	double p = getaveragep();
	ofs << ',' << p;
	getaveraget();
	double q = getaverageQ();
	ofs << ',' << q;
	double alfaK = q*alfaA + (1 - q)*alfaa;
	double betaK = q*betaA + (1 - q)*beta;
	double TPredicted = p*alfaK / (2 * betaK);
	ofs << ',' << TPredicted << std::endl;

}


void recordParameter() {
	para << "  Seed = " << seed << "\n";
	para << " iNgenesQ = " << ',' << iNgenesQ << "\n";
	para << " m = " << ',' << m << "\n";
	para << " ni = " << ',' << ni << "\n";
	para << " signma = " << ',' << signma << "\n";
	para << " alfaA = " << ',' << alfaA << "\n";
	para << " alfaa = " << ',' << alfaa << "\n";
	para << " viaA = " << ',' << viaA << "\n";
	para << " viaa = " << ',' << viaa << "\n";
	para << " betaA = " << ',' << betaA << "\n";
	para << " betaa = " << ',' << beta << "\n";
	para << " yeta = " << ',' << yeta << "\n";
	para << " PopSize = " << ',' << iPopsize << "\n";
	para << " Tmax = " << ',' << tmax << "\n";
	para << " f = " << ',' << f << "\n";
	para << " e = " << ',' << e << "\n";
	para << " emin = " << ',' << emin << "\n";
	para << " emax = " << ',' << emax << "\n";
	para << " Initial variation witdh = " << ',' << vartp << "\n";
	para << " Nsimulation = " << ',' << iNsimulation << "\n";
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

}

/*Main (se mai ci arriverò)(sono già qua) ( e siamo a 350 lines!)(450 now baby)*/


int main() {

	oo.open("Individual.csv", std::ios_base::app);
	ofs.open("Alltrait.csv", std::ios_base::app);
	para.open("Parameter.csv", std::ios_base::app);


	recordParameter();
	for (int ns = 0; ns < iNsimulation; ++ns) {
		// initial shit
		std::cout << " SIMULATION N =  " << ns << std::endl;
		oo << " \n Simulation N = " << ns << "\n";
		if (ns != 0)
			refresh();

		// initial shit

		oo << " \n Simulation N = " << seed << "\n";

		ofs << "  " << std::endl;
		ofs << " Seed " << ',' << " Time t " << ',' << " Preference p " << ',' << " Trait t1  " << ',' << " Trait t2 " << ',' << " t express " << ',' << " Mean quality" << ',' << " T predicted" << std::endl;


		// initial shit


		/*START SIMULATION FOR A LOT OF GENERATIONSSSSSS*/

		for (int t = 0; t < tmax; ++t) {


			if (t % 100 == 0) {
				oo << " \n NEXT GENERATION " << t << "\n";
				std::cout << " \n Generation n =  " << t << "\n";
				writeshit();
				ofs << seed << ',' << t;
				writeState();
			}



			// creat next generation
			NextGeneration();


			// Survivors became PARENTS!! congratss and CLEAR OFF POP, save p and t traits
			update();


		}
	}
	oo.close();
	ofs.close();
	para.close();

	return 0;
}


