/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package Bacterial_groups;
import java.io.File;
import svasta.FileRW;

/**
 *
 * @author Jelena
 */
public class MakeListOfOrganismsAndGenelements {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String path="C:\\Users\\Jelena\\Desktop\\transfer\\genbank\\genomes\\Bacteria";
        File file=new File(path);
        String path2="C:\\Users\\Jelena\\Desktop\\transfer\\ListOfOrganismsAndGenElementsAll.txt";
        FileRW fifi=new FileRW(path2, "W");
        if (file.isDirectory()) {
            File[] organismList = file.listFiles();
            for (int i = 0; i < organismList.length; i++) {
                File organism = organismList[i];

                //what is stored in an organism folder...
                if (organism.isDirectory()){
                    fifi.print(">"+organism.getName()+"\n", false);
                    //System.out.println("doing " + organismList[i]);
                    File[] fileList = organism.listFiles();
                    for (int k = 0; k < fileList.length; k++) {
                        //choosing to read in only .gbk files:
                        if(fileList[k].getName().contains(".gbk")){
                            FileRW f=new FileRW(fileList[k].getAbsolutePath(), "R");
                            String line;
                            String genelement=fileList[k].getName().substring(0, fileList[k].getName().length()-4);
                            String designation="";
                            String length="";
                            int z=0;
                            while((line=f.readLine())!=null){
                                if(line.startsWith("LOCUS")){
                                    String[] h=line.split("\\s+");
                                    length=h[2];
                                    z++;
                                }
                                else if(line.startsWith("DEFINITION")){
                                    String[] p=line.split(",");
                                    String[] h=p[0].split("\\s+");
                                    if(line.contains("genome")){
                                        designation="chromosome";
                                    }
                                    else if(line.contains("chromosome")){
                                        designation="chromosome"+h[h.length-1];
                                    }
                                    else if(line.contains("plasmid")){
                                        if(h[h.length-1].startsWith("p") || h[h.length-1].startsWith("mega")){
                                            designation=""+h[h.length-1];
                                        }
                                        else{
                                            designation="p"+h[h.length-1];
                                        }
                                    }
                                    else{
                                        System.out.println(line);
                                    }
                                    z++;
                                    
                                }
                                if(z==2){
                                    break;
                                }
                            }
                            fifi.print(genelement+"\t"+designation+"\t"+length+"\n", false);

                        }
                    }
                }
            }
            fifi.print("", true);
        }
        

    }

}

