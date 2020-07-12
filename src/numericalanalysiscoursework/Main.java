
package numericalanalysiscoursework;

/**
 *
 * @author Peter
 */

public class Main {

    public static double[][] classicalGS(double[][] inputVectors) {
        // Initialise
        double[][] outputVectors = new double[inputVectors.length][inputVectors[0].length];
        
        for (int i = 0; i < inputVectors.length; i++) {
            for (int j = 0; j < inputVectors[0].length; j++) {
                outputVectors[i][j] = inputVectors[i][j];
            }
        }
        
        for (int i = 0; i < inputVectors.length; i++) {
            for (int j = 0; j < i; j++) {
                double scalar = dotProduct(inputVectors[i], outputVectors[j]);
                        
                for (int k = 0; k < outputVectors[0].length; k++) {
                    outputVectors[i][k] -= scalar * outputVectors[j][k]; 
                }   
            }
            outputVectors[i] = normalise(outputVectors[i]);
        }
        
        return outputVectors;
    }
    
    public static double[][] modifiedGS(double[][] inputVectors) {
        // Initialise
        double[][] outputVectors = new double[inputVectors.length][inputVectors[0].length];
        
        // Initialise intermediate vectors
        double[][] interVectors = inputVectors;
        
        // Set q_1
        outputVectors[0] = normalise(interVectors[0]);
        
        for (int i = 1 ; i < inputVectors.length; i++) {
            for (int j = i; j < inputVectors.length; j++) {
                double scalar = dotProduct(interVectors[j], outputVectors[i - 1]);
                
                for (int k = 0; k < inputVectors[0].length; k++) {
                    // Change entries element by element
                    interVectors[j][k] -= (scalar * outputVectors[i - 1][k]);
                }
            }
            outputVectors[i] = normalise(interVectors[i]);
        }
        
        return outputVectors;
        
    }
        
    public static double dotProduct(double[] firstVector, double[] secondVector) {
        // Dot product
        double sum = 0;
        
        for (int i = 0; i < firstVector.length; i++) {
            sum += firstVector[i] * secondVector[i];
        }
        
        return sum;
    }
    
    public static double[] normalise(double[] vector) {
        // Calculate norm
        double norm = dotProduct(vector, vector);
        norm = 1 / Math.sqrt(norm);
        // Normalise
        vector = scalarMult(vector, norm);
        
        return vector;
    }
    
    public static double[] scalarMult(double[] vector, double scalar) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = vector[i] * scalar;
        }
        
        return vector;
    }
    
    public static void printMatrix(double[][] matrix) {
        // Method for printing a general matrix
        
        for (int i = 0; i < matrix[0].length; i++) {
            String str = "";
            for (int j = 0; j < matrix.length; j++) {
                str += matrix[j][i];
                
                if (j != matrix.length - 1) {
                    str += "\t";
                }
            }
            System.out.println(str);
        }
    }

    public static double[][] transpose(double[][] input) {
        // Method for transposing a matrix
        double[][] output = new double[input[0].length][input.length];
        
        for (int i = 0; i < input.length; i++) {
            for (int j = 0; j < input[0].length; j++) {
                output[j][i] = input[i][j];
            }
        }
        
        return output;
    }
    
    public static double[][] matrixMult(double[][] A, double[][] B) {
        
        A = transpose(A);
        B = transpose(B);

        int aRows = A.length;
        int aCols = A[0].length;
        int bRows = B.length;
        int bCols = B[0].length;

        double[][] C = new double[aRows][bCols];
        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bCols; j++) {
                C[i][j] = 0.0;
            }
        }
        
        // Abort if dimensions incorrect
        if (aCols != bRows) {
            throw new IllegalArgumentException("Dimensions not compatible");
        }

        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bCols; j++) {
                for (int k = 0; k < aCols; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        C = transpose(C);
        
        return C;
    }
        
    public static void main(String[] args) {
        
        // Convention: Matrices are written as columns: {{column one}, {column two}, ... }
        
        // Examples 

        double k = Math.pow(10, -10);
        double[][] inputVectors = new double[][]{
            {0, 0, 0, k, 1},
            {0, 0, k, 0, 1},
            {0, k, 0, 0, 1},
            {k, 0, 0, 0, 1}
        };
        
        double[][] identity = new double[][] {
            {1, 0, 0, 0, 0},
            {0, 1, 0, 0, 0},
            {0, 0, 1, 0, 0},
            {0, 0, 0, 1, 0},
            {0, 0, 0, 0, 1}
        };
        
        
        double[][] clasVectors = classicalGS(inputVectors);
        
        System.out.println("--- Q obtained using Classical Gram-Schmidt Algorithm ---");
        printMatrix(clasVectors);
        
        double[][] modfVectors = modifiedGS(inputVectors);
        
        System.out.println();
        System.out.println("--- Q bar obtained using Modified Gram-Schmidt Algorithm ---");
        printMatrix(modfVectors);
        
        System.out.println();
        System.out.println("--- Using CGS: Q^T Q - I^(n) ---");
        
        double[][] clasMatrix = matrixMult(transpose(clasVectors), clasVectors);
        
        for (int i = 0; i < clasMatrix.length; i++) {
            clasMatrix[i][i] -= identity[i][i];
        }
        
        printMatrix(clasMatrix);
        
        System.out.println();
        System.out.println("--- Using MGS: Q^T Q - I^(n) [where Q is Q bar] ---");
        
        double[][] modfMatrix = matrixMult(transpose(modfVectors), modfVectors);
        
        for (int i = 0; i < modfMatrix.length; i++) {
            modfMatrix[i][i] -= identity[i][i];
        }
        
        printMatrix(modfMatrix);
        
    }
    
}