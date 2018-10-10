package test;

import java.util.Vector;

/**
 * Producer-Consumer concurrency pattern. Based on code witten by Javin Paul:
 * <a href="http://www.java67.com/2012/12/producer-consumer-problem-with-wait-and-notify-example.html">
 * http://www.java67.com/2012/12/producer-consumer-problem-with-wait-and-notify-example.html</a> <br>
 * 
 * 
 * @author Dan Fulea, 18 AUG. 2016
 *
 */
public class ProducerConsumerDemo {
	public static void main(String args[]) {
        Vector<Integer> sharedQueue = new Vector<Integer>();
        int size = 1;//1;//4;
        int maxValue = 7;
        Thread prodThread = new Thread(new Producer(sharedQueue, size, maxValue), "Producer");
        Thread consThread = new Thread(new Consumer(sharedQueue, maxValue), "Consumer");
        prodThread.start();
        consThread.start();
    }
}
///////////ATOMIC==================
//Here's an example, because an example is often clearer than a long explanation.
//Suppose foo is a variable of type long. The following operation is not an atomic
//operation:

//foo = 65465498L;
//Indeed, the variable is written using two separate operations: one that 
//writes the first 32 bits, and a second one which writes the last 32 bits. 
//That means that another thread might read the value of foo, and see the intermediate 
//state. Making the operation atomic consists in using synchronization mechanisms in
//order to make sure that the operation is seen, from any other thread, as a single, 
//atomic (i.e. not splittable in parts), operation. That means that any other thread, 
//once the operation is made atomic, will either see the value of foo before the
//assignment, or after the assignment. But never the intermediate value.

//A simple way of doing this is to make the variable volatile:

//private volatile long foo;
//or to synchronize every access to the variable:

//public synchronized void setFoo(long value) {
//    this.foo = value;
//}

//public synchronized void getFoo() {
 //   return this.foo;
//}
// no other use of foo outside of these two methods, unless also synchronized
//Or to replace it with an AtomicLong:

//private AtomicLong foo;


//The java.util.Vector methods are all synchronized. So using it from multiple threads
//is "safe". You only need to synchronize if you need a read-evaluate-write process to be 
//atomic. Synchronizing your own methods does not necessarily make your code thread-safe
//for those scenarios. If the shared state is the Vector object, then you need to 
//synchronize on the instance of your Vector, not on an instance of your own classes.

//As stated above, every single method of Vector is thread-safe by its own because of 
//synchronized modifiers. But, if you need some complex operations, such as get() or 
//add() based on condition which is related to the same vector, this is not thread-safe.
//See example below:

//if (vector.size() > 0) {
 //   System.out.println(vector.get(0));
//}
//This code has a race condition between size() and get() - the size of vector might be 
//changed by other thread after our thread verified the vector is not empty, and 
//thus get() call may return unexpected results. To avoid this, the sample above 
//should be changed like this:

//synchronized (vector) {
//    if (vector.size() > 0) {
//        System.out.println(vector.get(0));
//    }
//}
//Now this "get-if-not-empty" operation is atomic and race condition-free.

//=========================================
class Producer implements Runnable {

    private final Vector<Integer> sharedQueue;
    private final int SIZE;
    private final int maxValue;

    public Producer(Vector<Integer> sharedQueue, int size, int maxValue) {
        this.sharedQueue = sharedQueue;
        this.SIZE = size;
        this.maxValue = maxValue;
    }

   
    public void run() {
        for (int i = 0; i < maxValue; i++) {
            //System.out.println("Produced: " + i);
            try {
                produce(i);
            } catch (InterruptedException ex) {
                ex.printStackTrace();
            }

        }
    }

    private void produce(int i) throws InterruptedException {

        //wait if queue is full
        while (sharedQueue.size() == SIZE) {
            synchronized (sharedQueue) {
                System.out.println("Queue is full " + Thread.currentThread().getName()
                                    + " is waiting , size: " + sharedQueue.size());

                sharedQueue.wait();//You need to be in a synchronized block in order for Object.wait() to work.
            }
        }

        //producing element and notify consumers
        synchronized (sharedQueue) {
            sharedQueue.add(i);
            System.out.println("Produced: " + i);
            sharedQueue.notifyAll();
        }
    }
}

class Consumer implements Runnable {

    private final Vector<Integer> sharedQueue;
    private final int maxValue;
    private boolean stopB=false;
    
    public Consumer(Vector<Integer> sharedQueue, int maxValue) {
        this.sharedQueue = sharedQueue;
        this.maxValue = maxValue;
    }

    
    public void run() {
        while (true && !stopB) {//!!!!!!!!!!!!!!!!
            try {
                //System.out.println("Consumed: " + consume());
            	consume();
                Thread.sleep(50);
            } catch (InterruptedException ex) {
                ex.printStackTrace();
            }

        }
    }

    private int consume() throws InterruptedException {
        //wait if queue is empty
        while (sharedQueue.isEmpty()){// && !stopB) {
            synchronized (sharedQueue) {
                System.out.println("Queue is empty " + Thread.currentThread().getName()
                                    + " is waiting , size: " + sharedQueue.size());

                sharedQueue.wait();
            }
        }

        //Otherwise consume element and notify waiting producer
        
        synchronized (sharedQueue) {
            sharedQueue.notifyAll();
            
            //==================================================
    	   
            if(sharedQueue.elementAt(0)==maxValue-1){
            	stopB = true;
            } else {
            	stopB=false;
            }
            //===================================================
           System.out.println("Consumed: " + sharedQueue.elementAt(0));
           return (Integer) sharedQueue.remove(0);
           
        }
    }
}
