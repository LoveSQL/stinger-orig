import java.util.Random;

public class main {
  static {
    System.loadLibrary("stinger");
  }

  public static void main(String[] Args) {
    int nv = 262144;
    int ne = nv * 8;

    int64Array sv = new int64Array(ne);
    int64Array ev = new int64Array(ne);
    int64Array w = new int64Array(ne);

    System.out.println("generating edges");
    Random rand = new Random();

    for(int i = 0; i < ne; i++) {
      sv.setitem(i, rand.nextInt(nv));
      ev.setitem(i, rand.nextInt(nv));
      w.setitem(i, 1);
    }

    System.out.println("creating stinger");
    SWIGTYPE_p_stinger s = stinger.edge_list_to_stinger(nv, ne, sv.cast(), ev.cast(), w.cast(), null, null, 0);

    System.out.println("new iterator");
    stinger_iterator_t it = stinger.stinger_iterator_new(s);
    int64Array types = new int64Array(1);
    types.setitem(0,0);

    System.out.println("setup iterator");
    stinger.stinger_iterator_edge_type_filter(types.cast(), 1, it);

    System.out.println("setting up components");
    boolean changed = false;
    int [] components = new int[nv];
    for(int i = 0; i < nv; i++) {
      components[i] = i;
    }

    System.out.println("calculating components");
    long start = System.currentTimeMillis();
    while(true) {
      changed = false;
      while(stinger.stinger_iterator_next(it) != 0) {
	if(components[(int)it.getSource()] < components[(int)it.getDest()]) {
	  components[(int)it.getDest()] = components[(int)it.getSource()];
	  changed = true;
	}
      }
      if(!changed) 
	break;
      for(int i = 0; i < nv; i++) {
	while(components[i] != components[components[i]]) {
	  components[i] = components[components[i]];
	}
      }
    }
    long stop = System.currentTimeMillis();

    System.out.println("done. time = " + ((double)(stop - start))/1000.0f);

    System.out.println("counting components");
    int count = 0;
    for(int i = 0; i < nv; i++) {
      if(components[i] == i)
	count++;
    }

    System.out.println(count);
  }
}
