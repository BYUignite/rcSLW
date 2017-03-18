//================================================================================
// This program uses the correlation and tabulated data to compute the ALBDF of a 
// single gas for given conditions.
//
// John Pearson
// September 20, 2013
//================================================================================

# include<iostream>
# include<math.h>
# include<cmath>
# include<fstream>
# include<vector>
# include<string>
using namespace std;

// Function Declarations
void data_read(vector<double>& Fdata1, vector<double>& Fdata2, int molecule, int P1);
double TABLE(double Cabs, double Tg, double Tb, double Y, double P, int P1, int molecule, vector<double>& C, vector<double>& Fdata1, vector<double>& Fdata2);
void locate(vector<double>& xx, int n, double x, int& j);
double ALBDF(double Cabs, double Tg, double Tb, double Y, double P, int molecule);

int main()
{

vector<double> Fdata1(2), Fdata2(2), Ctab(71), PPP(10);
double Tg, Tb, P, Y, C, ii;
double F_correlation, F_tables;
int P1, i, molecule;

/* ---------------------------------------------------------------------------- */
// Conditions at which to calculate the ALBDF

molecule=1; // Identify the species of interest: 1 for H2O, 2 for CO2, 3 for CO
C=1;		// m2/mole - the value of absorption cross-section
Tg=1500;	// K (gas temperature)
Tb=1500; 	// K (blackbody source temperature)
P=2;	 	// atm (total pressure)
Y=0.4; 		// mole fraction of the species of interest 
			// (only used for the case of H2O)

/* ---------------------------------------------------------------------------- */

// Set up absorption cross-section data for tabulated data
for (i=0; i<71; i++)
{
	ii=i;
	Ctab[i]=1e-4*pow((1000/1e-4),(ii/(70)));
}

// Get index of lower pressure bound
PPP[0]=0.1; PPP[1]=0.25; PPP[2]=0.5; PPP[3]=1; PPP[4]=2;
PPP[5]=4; PPP[6]=8; PPP[7]=15; PPP[8]=30; PPP[9]=50;
locate(PPP,9,P,P1);

// Read data files needed
data_read(Fdata1,Fdata2,molecule,P1);

// Calculate ALBDF using look-up table
F_tables=TABLE(C, Tg, Tb, Y, P, P1, molecule, Ctab, Fdata1, Fdata2);

// Calculate ALBDF using correlation
F_correlation=ALBDF(C, Tg, Tb, Y, P, molecule);

// Output results to file
ofstream fileout;

fileout.open("ALBDF.txt");

fileout << "F_table  F_correlation" << endl;
fileout << F_tables << "  " << F_correlation << endl;

fileout.close();

return(0);
}

//================================================================================
// This function computes the Absorption Line Blackbody Distribution Function
// (ALBDF) using the correlation.  Shown here is the case where total pressures
// other than atmospheric are expected.

double ALBDF(double Cabs, double Tg, double Tb, double Y, double PT, int molecule)
{
	double b[64], p[64], xi, xi_p, P, P2, F;
	int i, l, m, n;
	ifstream bfile, pfile;

	if (molecule==1) // H2O
	{	
		b[0]=0.978924033; b[1]=-0.207425641; b[2]=-2.806734007; b[3]=1.389417761;
		b[4]=2.43339282; b[5]=-5.542124572; b[6]=16.93794121; b[7]=-7.319199757; 
		b[8]=-2.270430289; b[9]=9.802027047; b[10]=-24.43728968; b[11]=10.2660654;
		b[12]=0.968416369; b[13]=-4.416640346; b[14]=10.26609154; b[15]=-4.203825171;
		b[16]=0.231177379; b[17]=-0.094605006; b[18]=0.908439931; b[19]=-0.538832358;
		b[20]=0.253813018; b[21]=2.082040478; b[22]=-4.476709604; b[23]=2.443789298; 
		b[24]=-0.337418966; b[25]=-2.308644559; b[26]=4.634251674; b[27]=-2.473496475; 
		b[28]=0.179162546; b[29]=0.78797212; b[30]=-1.520752749; b[31]=0.829104247;
		b[32]=0.166909096; b[33]=-0.795281465; b[34]=0.93926916; b[35]=-0.304702655;
		b[36]=-0.669061231; b[37]=3.939405065; b[38]=-4.430494571; b[39]=1.196316393;
		b[40]=0.83751338; b[41]=-4.931136293; b[42]=5.027576374; b[43]=-0.976626716;
		b[44]=-0.321006927; b[45]=1.907332712; b[46]=-1.76313509; b[47]=0.215566922;
		b[48]=0.018988364; b[49]=-0.075587929; b[50]=0.069280395; b[51]=2.48E-06;
		b[52]=-0.079441343; b[53]=0.346149251; b[54]=-0.231488972; b[55]=-0.065735266;
		b[56]=0.096930043; b[57]=-0.404603524; b[58]=0.151771109; b[59]=0.17764551;
		b[60]=-0.036760921; b[61]=0.147106808; b[62]=-0.014863816; b[63]=-0.098469241;

		p[0]=-15.032; p[1]=41.16; p[2]=-44.062; p[3]=15.958;
		p[4]=2.8319; p[5]=-1.4503; p[6]=4.1039; p[7]=-2.7305;
		p[8]=-0.17598; p[9]=2.0712; p[10]=-0.95861; p[11]=-0.16164;
		p[12]=-0.080753; p[13]=0.42606;	p[14]=-0.41713; p[15]=0.12426;
		p[16]=51.46; p[17]=-129.85;	p[18]=128.98; p[19]=-43.877;
		p[20]=-10.451; p[21]=11.362; p[22]=-31.935; p[23]=19.613;
		p[24]=2.4344; p[25]=-13.24; p[26]=7.4855; p[27]=0.030561;
		p[28]=0.56745; p[29]=-2.6102; p[30]=2.6208; p[31]=-0.807;
		p[32]=-49.884; p[33]=105.44; p[34]=-85.162; p[35]=22.81;
		p[36]=11.367; p[37]=-24.248; p[38]=66.901; p[39]=-39.492;
		p[40]=-5.9424; p[41]=24.356; p[42]=-14.925; p[43]=0.80233;
		p[44]=-1.0839; p[45]=4.5969; p[46]=-4.6497; p[47]=1.4447;
		p[48]=14.967; p[49]=-21.834; p[50]=6.4337; p[51]=2.4764;
		p[52]=-4.1549; p[53]=14.831; p[54]=-39.792; p[55]=22.768;
		p[56]=3.6337; p[57]=-13.111; p[58]=8.1435; p[59]=-0.5727;
		p[60]=0.57955; p[61]=-2.3468; p[62]=2.3494; p[63]=-0.72409;
	}

	if (molecule==2) // CO2
	{
		b[0]=1.834920165; b[1]=-1.822182621; b[2]=1.278328751; b[3]=-0.20979996;
		b[4]=-1.549692548; b[5]=1.11311754; b[6]=-2.425810746; b[7]=0.919085662;
		b[8]=3.517875983; b[9]=-1.149546992; b[10]=3.546566503; b[11]=-1.516418994;
		b[12]=-1.505796719; b[13]=0.522031728; b[14]=-1.653998778; b[15]=0.729684329;
		b[16]=0.202156834; b[17]=0.885000092; b[18]=-0.982474572; b[19]=0.375314671;
		b[20]=-0.531916205; b[21]=-2.214482516; b[22]=2.293219643; b[23]=-0.867071711;
		b[24]=1.019794785; b[25]=1.912899157; b[26]=-1.688329914; b[27]=0.600486161;
		b[28]=-0.506065247; b[29]=-0.52642387; b[30]=0.340073377; b[31]=-0.101266268;
		b[32]=0.054472216; b[33]=-0.1581808; b[34]=0.29003629; b[35]=-0.137608476;
		b[36]=-0.172336869; b[37]=0.538060901; b[38]=-0.853519053; b[39]=0.383655517;
		b[40]=0.233085884; b[41]=-0.702110304; b[42]=1.046918448; b[43]=-0.458832711;
		b[44]=-0.10143316; b[45]=0.29115141; b[46]=-0.424020704; b[47]=0.184047896;
		b[48]=0.006301556; b[49]=-0.030097135; b[50]=0.044577943; b[51]=-0.019483088;
		b[52]=-0.018227338; b[53]=0.095092466; b[54]=-0.122626495; b[55]=0.05015993;
		b[56]=0.020997129; b[57]=-0.109841221; b[58]=0.134896941; b[59]=-0.05405009;
		b[60]=-8.18E-03; b[61]=0.04163713; b[62]=-0.05021366; b[63]=0.020023722;

		p[0]=-14.9838124; p[1]=38.41682512; p[2]=-37.08708237; p[3]=11.6339692;
		p[4]=-0.279676701; p[5]=-2.38648754; p[6]=3.488321976; p[7]=-1.655314601;
		p[8]=0.240326019; p[9]=-1.253645566; p[10]=0.926603383; p[11]=-0.175741918;
		p[12]=0.026838401; p[13]=-0.074737215; p[14]=0.010782536; p[15]=0.020039736;
		p[16]=65.6646184; p[17]=-153.1300705; p[18]=134.5182112; p[19]=-38.49884117;
		p[20]=3.731942483; p[21]=-6.345566639; p[22]=4.113471062; p[23]=0.442687983;
		p[24]=-1.457123086; p[25]=4.493581565; p[26]=-1.904072879; p[27]=-0.230432069;
		p[28]=-0.194708661; p[29]=0.344715216; p[30]=0.082301336; p[31]=-0.169291984;
		p[32]=-90.49778524; p[33]=191.2774078; p[34]=-145.3015474; p[35]=34.01824433;
		p[36]=-10.43915684; p[37]=40.10843718; p[38]=-43.29133973; p[39]=12.82988374;
		p[40]=2.551054114; p[41]=-6.770847243; p[42]=2.075279528; p[43]=0.815010456;
		p[44]=0.400183431; p[45]=-0.766811133; p[46]=0.018447262; p[47]=0.238835363;
		p[48]=41.94117422; p[49]=-79.99441926; p[50]=49.60377935; p[51]=-7.294292574;
		p[52]=7.51406621; p[53]=-36.3335072; p[54]=42.51857862; p[55]=-14.30274633;
		p[56]=-1.431539415; p[57]=4.097348243; p[58]=-1.864356794; p[59]=-0.109473928;
		p[60]=-0.256915604; p[61]=0.610041193; p[62]=-0.249699442; p[63]=-0.038328249;
	}

	if (molecule==3) // CO
	{
		b[0]=3.316310112; b[1]=0.273254392; b[2]=-1.035375039; b[3]=0.473979958;
		b[4]=-7.974114073; b[5]=0.546750056; b[6]=1.639880389; b[7]=-0.875343857;
		b[8]=15.32452431; b[9]=-3.644252282; b[10]=-0.372777673; b[11]=0.687532486;
		b[12]=-7.962324128; b[13]=3.982035927; b[14]=-2.024788009; b[15]=0.457744088;
		b[16]=0.427224963; b[17]=0.159347829; b[18]=-0.314631881; b[19]=0.152076858;
		b[20]=-1.356291059; b[21]=0.115369369; b[22]=0.658528303; b[23]=-0.372592736;
		b[24]=3.663861248; b[25]=-1.945483178; b[26]=0.377631485; b[27]=0.088835684;
		b[28]=-2.530222444; b[29]=2.519532063; b[30]=-1.79410425; b[31]=0.564326917;
		b[32]=0.072745637; b[33]=-0.061526669; b[34]=0.033240467; b[35]=-0.004119976;
		b[36]=-0.266116052; b[37]=0.003566603; b[38]=0.127097068; b[39]=-0.066171379;
		b[40]=0.704482493; b[41]=-0.253535518; b[42]=-0.009919421; b[43]=0.031869833;
		b[44]=-0.491444366; b[45]=0.432780135; b[46]=-0.32227219; b[47]=0.111101186;
		b[48]=0.004004596; b[49]=-0.005870558; b[50]=0.004038598; b[51]=-0.000880491;
		b[52]=-0.015436392; b[53]=-0.00290339; b[54]=0.010857925; b[55]=-0.005012608;
		b[56]=0.042245518; b[57]=-0.007520854; b[58]=-0.007585042; b[59]=0.003917497;
		b[60]=-0.029663594; b[61]=0.022352969; b[62]=-0.016590884; b[63]=0.006069374;

		p[0]=-7.980270361; p[1]=12.44840073; p[2]=-13.11253512; p[3]=4.868721547;
		p[4]=-1.029852518; p[5]=5.343062546; p[6]=-4.232122745; p[7]=1.063879591;
		p[8]=-0.12552327; p[9]=0.29017319; p[10]=0.608824895; p[11]=-0.457268847;
		p[12]=-0.010648303; p[13]=-0.008657971; p[14]=0.096011021; p[15]=-0.056109023;
		p[16]=30.9289566; p[17]=-62.58487087; p[18]=76.96827504; p[19]=-30.93824358;
		p[20]=7.656923091; p[21]=-28.79840669; p[22]=26.67286357; p[23]=-8.605697097;
		p[24]=2.104107011; p[25]=-2.63598025; p[26]=-1.777881616; p[27]=1.651348405;
		p[28]=0.177021253; p[29]=-0.05593256; p[30]=-0.401553289; p[31]=0.245224992;
		p[32]=-34.73465418; p[33]=103.8477015; p[34]=-141.729319; p[35]=58.45442338;
		p[36]=-15.4237926; p[37]=49.47310879; p[38]=-48.57742671; p[39]=16.9208218;
		p[40]=-5.547243446; p[41]=5.732398138; p[42]=1.987247698; p[43]=-2.23242013;
		p[44]=-0.466208091; p[45]=0.190251348; p[46]=0.635513731; p[47]=-0.389171949;
		p[48]=11.61823155; p[49]=-57.78495547; p[50]=84.50680302; p[51]=-35.10802171;
		p[52]=8.766127928; p[53]=-27.37822507; p[54]=28.39218696; p[55]=-10.34802912;
		p[56]=3.655832622; p[57]=-3.373560674; p[58]=-0.873746551; p[59]=1.050506404;
		p[60]=0.305051681; p[61]=-0.094529018; p[62]=-0.378666546; p[63]=0.219049669;
	}

	// Calculate pressure term
	xi=log(Cabs);
	xi_p=0;

	// equivalent pressure for H2O
	if (molecule==1)
		PT=(1+8.17*Y)*PT;

	i=0;
	for (l=0; l<4; l++) 
	{
		for (m=0; m<4; m++) 
		{
			for (n=0; n<4; n++) 
			{

				xi_p=xi_p+p[i]*pow(Tg*Tb/(2500*2500),n)*pow(xi,m)*pow(log(PT*100)/10,l+1);
				i++;
			}
		}
	}

	// Calculate P
	P=0;
	P2=0;
	
	if (molecule==3)
	{
		P2=xi_p;
		xi_p=0;
	}

	i=0;
	for (l=0; l<4; l++)
	{
		for (m=0; m<4; m++)
		{
			for (n=0; n<4; n++)
			{
			
				P=P+b[i]*pow((Tg/2500),n)*pow((Tb/2500),m)*pow((xi-xi_p),l);
				i++;
			}
		}
	}

	// Calculate ALBDF
	F=.5*tanh(P-P2)+.5;

	return F;
}


//================================================================================
// This function computes the Absorption Line Blackbody Distribution Function 
// (ALBDF) using tabulated data.

double TABLE(double Cabs, double Tg, double Tb, double Y, double P, int P1, int molecule, vector<double>& C, vector<double>& Fdata1, vector<double>& Fdata2)
{

int i, j, k, l, m;
double F;
vector<double> T(28), Fint(16), YY(9), PP(10);

// Put values in range if they are out of bounds because
// extrapolation is not allowed when using tabulated data
if (Cabs<C[0]) Cabs=C[0];
if (Cabs>C[70]) Cabs=C[70];
if (Tg<300) Tg=300;
if (Tg>3000) Tg=3000;
if (Tb<300) Tb=300;
if (Tb>3000) Tb=3000;

// Set temperature values
T[0]=300;
for (i=1; i<28; i++)
	T[i]=T[i-1]+100;

// Set pressure values
PP[0]=0.1; PP[1]=0.25; PP[2]=0.5; PP[3]=1; PP[4]=2;
PP[5]=4; PP[6]=8; PP[7]=15; PP[8]=30; PP[9]=50;
	
// Find lower indices needed for interpolation in Tg, Tb, and C
locate(T, 27, Tg, m);
locate(T, 27, Tb, l);
locate(C, 70, Cabs, k);

// H2O section
if (molecule==1)
{

	// Set mole fraction values
	YY[0]=0; YY[1]=0.05; YY[2]=0.1; YY[3]=0.2; YY[4]=0.3;
	YY[5]=0.4; YY[6]=0.6; YY[7]=0.8; YY[8]=1;

	// Find lower index needed for interpolation in Y
	locate(YY, 8, Y, j);
	
	// Beginning index for interpolation
	i=j*55664+m*1988+l*71+k;
	
	// Another way to perform this interpolation would be to use a multi-dimensional vector
	// for the tabulated data instead of a single column i.e. Fdata[N_C][N_Tg][N_Tb][N_Y][N_P].
	// For this method, setting up the interpolation would be simpler. However, locating 
	// the needed points for the interpolation when using a single column of data is also 
	// straightforward and is demonstrated here.

	// Interpolate in P - linear interpolation applied
	Fint[0]=Fdata1[i]+(Fdata2[i]-Fdata1[i])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[1]=Fdata1[i+55664]+(Fdata2[i+55664]-Fdata1[i+55664])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[2]=Fdata1[i+1]+(Fdata2[i+1]-Fdata1[i+1])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[3]=Fdata1[i+55665]+(Fdata2[i+55665]-Fdata1[i+55665])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[4]=Fdata1[i+71]+(Fdata2[i+71]-Fdata1[i+71])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[5]=Fdata1[i+55735]+(Fdata2[i+55735]-Fdata1[i+55735])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[6]=Fdata1[i+72]+(Fdata2[i+72]-Fdata1[i+72])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[7]=Fdata1[i+55736]+(Fdata2[i+55736]-Fdata1[i+55736])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[8]=Fdata1[i+1988]+(Fdata2[i+1988]-Fdata1[i+1988])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[9]=Fdata1[i+57652]+(Fdata2[i+57652]-Fdata1[i+57652])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[10]=Fdata1[i+1989]+(Fdata2[i+1989]-Fdata1[i+1989])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[11]=Fdata1[i+57653]+(Fdata2[i+57653]-Fdata1[i+57653])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[12]=Fdata1[i+2059]+(Fdata2[i+2059]-Fdata1[i+2059])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[13]=Fdata1[i+57723]+(Fdata2[i+57723]-Fdata1[i+57723])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[14]=Fdata1[i+2060]+(Fdata2[i+2060]-Fdata1[i+2060])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[15]=Fdata1[i+57724]+(Fdata2[i+57724]-Fdata1[i+57724])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	
	// Interpolate in Y
	Fint[0]=Fint[0]+(Fint[1]-Fint[0])*(Y-YY[j])/(YY[j+1]-YY[j]);
	Fint[1]=Fint[2]+(Fint[3]-Fint[2])*(Y-YY[j])/(YY[j+1]-YY[j]);
	Fint[2]=Fint[4]+(Fint[5]-Fint[4])*(Y-YY[j])/(YY[j+1]-YY[j]);
	Fint[3]=Fint[6]+(Fint[7]-Fint[6])*(Y-YY[j])/(YY[j+1]-YY[j]);
	Fint[4]=Fint[8]+(Fint[9]-Fint[8])*(Y-YY[j])/(YY[j+1]-YY[j]);
	Fint[5]=Fint[10]+(Fint[11]-Fint[10])*(Y-YY[j])/(YY[j+1]-YY[j]);
	Fint[6]=Fint[12]+(Fint[13]-Fint[12])*(Y-YY[j])/(YY[j+1]-YY[j]);
	Fint[7]=Fint[14]+(Fint[15]-Fint[14])*(Y-YY[j])/(YY[j+1]-YY[j]);

	// Interpolate in C
	Fint[0]=Fint[0]+(Fint[1]-Fint[0])*(Cabs-C[k])/(C[k+1]-C[k]);
	Fint[1]=Fint[2]+(Fint[3]-Fint[2])*(Cabs-C[k])/(C[k+1]-C[k]);
	Fint[2]=Fint[4]+(Fint[5]-Fint[4])*(Cabs-C[k])/(C[k+1]-C[k]);
	Fint[3]=Fint[6]+(Fint[7]-Fint[6])*(Cabs-C[k])/(C[k+1]-C[k]);

	// Interpolate in Tb
	Fint[0]=Fint[0]+(Fint[1]-Fint[0])*(Tb-T[l])/(T[l+1]-T[l]);
	Fint[1]=Fint[2]+(Fint[3]-Fint[2])*(Tb-T[l])/(T[l+1]-T[l]);

	// Interpolate in Tg
	F=Fint[0]+(Fint[1]-Fint[0])*(Tg-T[m])/(T[m+1]-T[m]);
}

// CO2 or CO section
if (molecule>1)
{

	i=m*1988+l*71+k;

	// Interpolate in P
	Fint[0]=Fdata1[i]+(Fdata2[i]-Fdata1[i])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[1]=Fdata1[i+1]+(Fdata2[i+1]-Fdata1[i+1])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[2]=Fdata1[i+71]+(Fdata2[i+71]-Fdata1[i+71])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[3]=Fdata1[i+72]+(Fdata2[i+72]-Fdata1[i+72])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[4]=Fdata1[i+1988]+(Fdata2[i+1988]-Fdata1[i+1988])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[5]=Fdata1[i+1989]+(Fdata2[i+1989]-Fdata1[i+1989])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[6]=Fdata1[i+2059]+(Fdata2[i+2059]-Fdata1[i+2059])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	Fint[7]=Fdata1[i+2060]+(Fdata2[i+2060]-Fdata1[i+2060])*(P-PP[P1])/(PP[P1+1]-PP[P1]);
	
	// Interpolate in C
	Fint[0]=Fint[0]+(Fint[1]-Fint[0])*(Cabs-C[k])/(C[k+1]-C[k]);
	Fint[1]=Fint[2]+(Fint[3]-Fint[2])*(Cabs-C[k])/(C[k+1]-C[k]);
	Fint[2]=Fint[4]+(Fint[5]-Fint[4])*(Cabs-C[k])/(C[k+1]-C[k]);
	Fint[3]=Fint[6]+(Fint[7]-Fint[6])*(Cabs-C[k])/(C[k+1]-C[k]);

	// Interpolate in Tb
	Fint[0]=Fint[0]+(Fint[1]-Fint[0])*(Tb-T[l])/(T[l+1]-T[l]);
	Fint[1]=Fint[2]+(Fint[3]-Fint[2])*(Tb-T[l])/(T[l+1]-T[l]);

	// Interpolate in Tg
	F=Fint[0]+(Fint[1]-Fint[0])*(Tg-T[m])/(T[m+1]-T[m]);
}

return F;
}

//================================================================================
// This function reads in the tabulated ALBDF data for H2O, CO2, and CO

void data_read(vector<double>& Fdata1, vector<double>& Fdata2, int molecule, int P1)
{

ifstream data_F1, data_F2;

if (molecule==1)
{
	
	// Specify file locations and names (make sure file location is set correctly)
	if (P1==0) {
		data_F1.open("ALBDF_Tables\\h2o_p0_1.txt");
		data_F2.open("ALBDF_Tables\\h2o_p0_25.txt");
	}
	if (P1==1) {
		data_F1.open("ALBDF_Tables\\h2o_p0_25.txt");
		data_F2.open("ALBDF_Tables\\h2o_p0_5.txt");
	}
	if (P1==2) {
		data_F1.open("ALBDF_Tables\\h2o_p0_5.txt");
		data_F2.open("ALBDF_Tables\\h2o_p1.txt");
	}
	if (P1==3) {
		data_F1.open("ALBDF_Tables\\h2o_p1.txt");
		data_F2.open("ALBDF_Tables\\h2o_p2.txt");
	}
	if (P1==4) {
		data_F1.open("ALBDF_Tables\\h2o_p2.txt");
		data_F2.open("ALBDF_Tables\\h2o_p4.txt");
	}
	if (P1==5) {
		data_F1.open("ALBDF_Tables\\h2o_p4.txt");
		data_F2.open("ALBDF_Tables\\h2o_p8.txt");
	}
	if (P1==6) {
		data_F1.open("ALBDF_Tables\\h2o_p8.txt");
		data_F2.open("ALBDF_Tables\\h2o_p15.txt");
	}
	if (P1==7) {
		data_F1.open("ALBDF_Tables\\h2o_p15.txt");
		data_F2.open("ALBDF_Tables\\h2o_p30.txt");
	}
	if (P1>7) {
		data_F1.open("ALBDF_Tables\\h2o_p30.txt");
		data_F2.open("ALBDF_Tables\\h2o_p50.txt");
	}
	
	Fdata1.resize(500976);	
	Fdata2.resize(500976);	
	for (int i=0; i<500976; i++) {
		data_F1 >> Fdata1[i];
		data_F2 >> Fdata2[i];
	}
}

if (molecule==2)
{

	if (P1==0) {
		data_F1.open("ALBDF_Tables\\co2_p0_1.txt");
		data_F2.open("ALBDF_Tables\\co2_p0_25.txt");
	}
	if (P1==1) {
		data_F1.open("ALBDF_Tables\\co2_p0_25.txt");
		data_F2.open("ALBDF_Tables\\co2_p0_5.txt");
	}
	if (P1==2) {
		data_F1.open("ALBDF_Tables\\co2_p0_5.txt");
		data_F2.open("ALBDF_Tables\\co2_p1.txt");
	}
	if (P1==3) {
		data_F1.open("ALBDF_Tables\\co2_p1.txt");
		data_F2.open("ALBDF_Tables\\co2_p2.txt");
	}
	if (P1==4) {
		data_F1.open("ALBDF_Tables\\co2_p2.txt");
		data_F2.open("ALBDF_Tables\\co2_p4.txt");
	}
	if (P1==5) {
		data_F1.open("ALBDF_Tables\\co2_p4.txt");
		data_F2.open("ALBDF_Tables\\co2_p8.txt");
	}
	if (P1==6) {
		data_F1.open("ALBDF_Tables\\co2_p8.txt");
		data_F2.open("ALBDF_Tables\\co2_p15.txt");
	}
	if (P1==7) {
		data_F1.open("ALBDF_Tables\\co2_p15.txt");
		data_F2.open("ALBDF_Tables\\co2_p30.txt");
	}
	if (P1>7) {
		data_F1.open("ALBDF_Tables\\co2_p30.txt");
		data_F2.open("ALBDF_Tables\\co2_p50.txt");
	}
	
	Fdata1.resize(55664);	
	Fdata2.resize(55664);	
	for (int i=0; i<55664; i++) {
		data_F1 >> Fdata1[i];
		data_F2 >> Fdata2[i];
	}
}

if (molecule==3)
{

	if (P1==0) {
		data_F1.open("ALBDF_Tables\\co_p0_1.txt");
		data_F2.open("ALBDF_Tables\\co_p0_25.txt");
	}
	if (P1==1) {
		data_F1.open("ALBDF_Tables\\co_p0_25.txt");
		data_F2.open("ALBDF_Tables\\co_p0_5.txt");
	}
	if (P1==2) {
		data_F1.open("ALBDF_Tables\\co_p0_5.txt");
		data_F2.open("ALBDF_Tables\\co_p1.txt");
	}
	if (P1==3) {
		data_F1.open("ALBDF_Tables\\co_p1.txt");
		data_F2.open("ALBDF_Tables\\co_p2.txt");
	}
	if (P1==4) {
		data_F1.open("ALBDF_Tables\\co_p2.txt");
		data_F2.open("ALBDF_Tables\\co_p4.txt");
	}
	if (P1==5) {
		data_F1.open("ALBDF_Tables\\co_p4.txt");
		data_F2.open("ALBDF_Tables\\co_p8.txt");
	}
	if (P1==6) {
		data_F1.open("ALBDF_Tables\\co_p8.txt");
		data_F2.open("ALBDF_Tables\\co_p15.txt");
	}
	if (P1==7) {
		data_F1.open("ALBDF_Tables\\co_p15.txt");
		data_F2.open("ALBDF_Tables\\co_p30.txt");
	}
	if (P1>7) {
		data_F1.open("ALBDF_Tables\\co_p30.txt");
		data_F2.open("ALBDF_Tables\\co_p50.txt");
	}
	
	Fdata1.resize(55664);	
	Fdata2.resize(55664);	
	for (int i=0; i<55664; i++) {
		data_F1 >> Fdata1[i];
		data_F2 >> Fdata2[i];
	}
}

data_F1.close();
data_F2.close();

}

//================================================================================
// This function searches data using bisection. The code is
// an adaptation of a code from Numerical Recipes in C.
// xx is vector to search, n is max index, x is value to search 
// for, and j is lower index of the two points x is between.

void locate(vector<double>& xx, int n, double x, int& j)
{

unsigned long ju, jm, jl;
int ascnd;

jl=0;
ju=n+1;
ascnd=(xx[n] >= xx[1]);
while (ju-jl>1) {
	jm=(ju+jl) >> 1;
	if (x >= xx[jm] == ascnd)
		jl=jm;
	else
		ju=jm;
}
if (x == xx[0]) j=0;
else if (x==xx[n]) j=n-1;
else j=jl;

}
