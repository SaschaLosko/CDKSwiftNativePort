import java.io.*;
import java.util.*;
import org.openscience.cdk.*;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.io.*;
import org.openscience.cdk.isomorphism.matchers.*;
public class RGroupParity {
 public static void main(String[] args) throws Exception {
  System.out.println("{"); boolean firstFile=true;
  for(int n : new int[]{1,2,3,4,5,7}) {
   try(var reader=new RGroupQueryReader(new FileInputStream(args[0]+"/rgfile."+n+".mol"))) {
    var query=reader.read(new RGroupQuery(DefaultChemObjectBuilder.getInstance()));
    var histogram=new TreeMap<String,Integer>();
    for(var molecule:query.getAllConfigurations()) {
     var atoms=new ArrayList<String>(); var bonds=new ArrayList<String>();
     for(var a:molecule.atoms()) atoms.add(a.getSymbol()+":"+a.getFormalCharge()+":"+Objects.toString(a.getProperty(CDKConstants.REST_H),"nil"));
     for(var b:molecule.bonds()) {String a=b.getBegin().getSymbol(),z=b.getEnd().getSymbol(); if(a.compareTo(z)>0){String t=a;a=z;z=t;} bonds.add(a+":"+b.getOrder().numeric()+":"+z);}
     Collections.sort(atoms); Collections.sort(bonds);
     String key=String.join(",",atoms)+"|"+String.join(",",bonds);
     histogram.merge(key,1,Integer::sum);
    }
    if(!firstFile)System.out.println(",");firstFile=false;
    System.out.print("\""+n+"\":{");boolean first=true;
    for(var e:histogram.entrySet()){if(!first)System.out.print(",");first=false;System.out.print("\""+e.getKey()+"\":"+e.getValue());}
    System.out.print("}");
   }
  }
  System.out.println("\n}");
 }
}
