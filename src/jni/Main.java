public class Main {
    public static void main(String[] args) {
        // 1️⃣ Initialize GetelecInterface
        GetelecInterface getelec = new GetelecInterface("getelec.cfg", "default");

        // 2️⃣ Generate xFN values (equivalent to NumPy linspace(0.1, 0.5, 512))
        double[] xFN = new double[64];
        for (int i = 0; i < xFN.length; i++) {
            xFN[i] = 0.1 + i * (0.5 - 0.1) / (64 - 1);
        }

        // 3️⃣ Compute 1 / xFN and set the field
        double[] invXFN = new double[xFN.length];
        for (int i = 0; i < xFN.length; i++) {
            invXFN[i] = 1.0 / xFN[i];
        }
        getelec.setField(invXFN);

        // 4️⃣ Run the computation
        getelec.run(false);

        // 5️⃣ Retrieve and print current densities
        double[] currentDensities = getelec.getCurrentDensities();

        System.out.println("Current Densities:");
        for (double density : currentDensities) {
            System.out.println(density);
        }
    }
}