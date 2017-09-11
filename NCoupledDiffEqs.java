// The purpose of this class is to solve N coupled differential equations using the RK4 method.
// Start Date: 09/22/2014

public class NCoupledDiffEqs {

	public static int N = 4; // number of differential equations, this number is hardcoded so the program can be run at once. Should be changed to parameter in an OO implemenation.
	public static double[] y = new double[N]; // y: (y1,y2,...,yn)
	public static double[] diffEq = new double[N];  // this array will contain the coupled differential equations
	public static double deltaX = 0.001;    // integration step size
	public static double y2 = 1;
	public static double z = 1;
	public static double x = 0;
	public static int time = 0;

	public static void main(String[] args) {
		
		initialize();
		for (int i = 0; i < 100; i++) {
			RKSolveN();
		}
	}
	// initial conditions for solutions
	public static void initialize()
	{
		for(int i=0;i<N;i++)
		{
			y[i] = 1;
		}
	}
	
	public static double[] update = new double[N];
	
	public static double dydx(double x, int i) 
	{
		dydxUpdate(update);
		return diffEq[i];
	}

	// define differential equations here
	public static void dydxUpdate(double[] y) 
	{	
		diffEq[0] = -y[1]-y[0]*y[0];
		diffEq[1] = 2*y[0]-y[1];
		diffEq[2] = -y[2]+y[0]*y[1];
		diffEq[3] = -y[3] +y[0];
	}

	public static void RKSolveN() {
		int steps = 100;
		double x;

		double[] k1 = new double[N];
		double[] k2 = new double[N];
		double[] k3 = new double[N];
		double[] k4 = new double[N];
		
		for(int i=0;i<y.length;i++)
		{
			update[i] = y[i];
		}
		
		for (int i = 0; i < steps; i++) {
			x = i * deltaX;

			for (int j = 0; j < N; j++) {
				k1[j] = deltaX*dydx(x,j);
			}
			
			for (int j = 0; j < N; j++) {
				update[j] = y[j]+k1[j]/2.0;
			}
			
			for (int j = 0; j < N; j++) {
				k2[j] = deltaX*dydx(x+deltaX/2.0,j);
			}
			
			for (int j = 0; j < N; j++) {
				update[j] = y[j]+k2[j]/2.0;
			}
			
			for (int j = 0; j < N; j++) {
				k3[j] = deltaX*dydx(x+deltaX/2.0,j);
			}
			
			for (int j = 0; j < N; j++) {
				update[j] = y[j]+k3[j];
			}
			
			for (int j = 0; j < N; j++) {
				k4[j] = deltaX*dydx(x+deltaX,j);
			}
			

			
			for(int j=0;j<N;j++)
			{
				y[j] += (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])/6.0;
			}
			
			if(time==0)
			{
				System.out.print(time+" ");
				for(int j=0;j<N;j++)
				{	
					System.out.print(" "+y[j]);
				}
				System.out.println();
			}
			time++;
		}
		
		System.out.print(time+" ");
		for(int j=0;j<N;j++)
		{	
			System.out.print(" "+y[j]);
		}
		System.out.println();

	}

	public static double F(double x, double y, double z) {
		return -z - y * y;
	}

	public static double G(double x, double y, double z) {
		return 2 * y - z;
	}

	public static void RKSolve() {
		int steps = 100;
		double x;
		double k1, k2, k3, k4, l1, l2, l3, l4;

		for (int i = 0; i < steps; i++) {
			x = i * deltaX;
			// for loop goes here.
			k1 = deltaX * F(x, y2, z);
			l1 = deltaX * G(x, y2, z);
			k2 = deltaX * F(x + deltaX / 2.0, y2 + k1 / 2.0, z + l1 / 2.0);
			l2 = deltaX * G(x + deltaX / 2.0, y2 + k1 / 2.0, z + l1 / 2.0);
			k3 = deltaX * F(x + deltaX / 2.0, y2 + k2 / 2.0, z + l2 / 2.0);
			l3 = deltaX * G(x + deltaX / 2.0, y2 + k2 / 2.0, z + l2 / 2.0);
			k4 = deltaX * F(x + deltaX, y2 + k3, z + l3);
			l4 = deltaX * G(x + deltaX, y2 + k3, z + l3);

			if (time == 0){System.out.println(time + " " + y2 + " " + z);}

			y2 += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
			z += (l1 + 2 * l2 + 2 * l3 + l4) / 6.0;
			
		
			time++;
		}
		System.out.println(time + " " + y2 + " " + z);
	}

}
