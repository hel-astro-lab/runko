#from mpi4py.futures import MPIPoolExecutor
from mpi4py import MPI

import argparse
from parser import parse_input


from merge_nodes import merge_field_nodes
from merge_nodes import merge_analysis_nodes



# based on https://github.com/jbornschein/mpi4py-examples/blob/master/09-task-pull.py


#x0, x1, w = -2.0, +2.0, 640*2
#y0, y1, h = -1.5, +1.5, 480*2
#dx = (x1 - x0) / w
#dy = (y1 - y0) / h
#
#c = complex(0, 0.65)
#
#def julia(x, y):
#    z = complex(x, y)
#    n = 255
#    while abs(z) < 3 and n > 1:
#        z = z**2 + c
#        n -= 1
#    return n
#
#def julia_line(k):
#    line = bytearray(w)
#    y = y1 - k * dy
#    for j in range(w):
#        x = x0 + j * dx
#        line[j] = julia(x, y)
#    return line


def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

# Define MPI message tags
tags = enum('READY', 'DONE', 'EXIT', 'START')






class Tasker:

    def get_args(self):
        conf, fdir, args = parse_input()

        return conf, fdir, args


    def __init__(self):

        # Initializations and preliminaries
        self.comm = MPI.COMM_WORLD   # get MPI communicator object
        self.size = self.comm.size   # total number of processes
        self.rank = self.comm.rank   # rank of this process
        self.status = MPI.Status()   # get MPI status object

        # parse who is master
        self.master = False
        if self.rank == 0:
            self.master = True

        self.conf, self.fdir, args = self.get_args()
        if self.master:
            print("Tasker operating on dir:{} with {} workers".format(self.fdir, self.size))







    ##################################################
    # task lists to be executed
    lap_tasks = [ 'test_member', 'test_other_member' ]

    sim_tasks = [ 'test_global' ]



    ##################################################
    # internal variables

    work_to_do = True # distribute tasks until this is false
    lap = 0  # keep track of lap we are processing
    task_index = 0 #keep track of lap task we are processing


    # task feeder
    def get_new_task(self):
        
        # cycle among laps
        if self.task_index >= len(self.lap_tasks):
            self.lap += self.conf.interval
            self.task_index = 0

        task = {
                'name': self.lap_tasks[self.task_index],
                'args': self.lap,
                }
        self.task_index += 1


        if self.lap >= 1000:
            self.work_to_do = False

        return task



    # task scheduling; master process executes this
    def schedule_tasks(self):

        #task_index = 0
        num_workers = self.size - 1
        closed_workers = 0

        print("Master starting with %d workers" % num_workers)
        while closed_workers < num_workers:
            data = self.comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=self.status)
            source = self.status.Get_source()
            tag = self.status.Get_tag()
            if tag == tags.READY:

                # Worker is ready, so send it a task
                #if task_index < len(tasks):
                if self.work_to_do:

                    task = self.get_new_task()
                    
                    self.comm.send(task, dest=source, tag=tags.START)
                    print("Sending task {} to worker {}".format(task['name'], source))
                else:
                    self.comm.send(None, dest=source, tag=tags.EXIT)


            elif tag == tags.DONE:
                results = data

                print("Got data from worker %d" % source)

            elif tag == tags.EXIT:
                print("Worker %d exited." % source)
                closed_workers += 1

        print("Master finishing...")



    # worker that takes a task and works on it
    def work_on_tasks(self):

        # Worker processes execute code below
        name = MPI.Get_processor_name()
        print("I am a worker with rank %d on %s." % (self.rank, name))

        while True:
            self.comm.send(None, dest=0, tag=tags.READY)
            task = self.comm.recv(source=0, tag=MPI.ANY_TAG, status=self.status)
            tag = self.status.Get_tag()
            if tag == tags.START:

                # Do the work here
                method_to_call = getattr(self, task['name'])
                result = method_to_call( task['args'] )

                self.comm.send(result, dest=0, tag=tags.DONE)
            elif tag == tags.EXIT:
                break

        self.comm.send(None, dest=0, tag=tags.EXIT)


    def test_member(self, args):
        print("running test_member with args", args)

        return True

    def test_other_member(self, args):
        print("running test_other_member with args", args)

        return True

    def test_global(self, args):
        print("running test_global with args", args)

        return True



    def run(self):

        if self.master:
            self.schedule_tasks()
        else:
            self.work_on_tasks()




if __name__ == '__main__':

    tasker = Tasker()
    tasker.run()




    #with MPIPoolExecutor() as executor:
    #    image = executor.map(julia_line, range(h))
    #    with open('julia.pgm', 'wb') as f:
    #        f.write(b'P5 %d %d %d\n' % (w, h, 255))
    #        for line in image:
    #            f.write(line)

    #executor = MPIPoolExecutor()
    #future = executor.submit(pow, 321, 1234)
    #print(future.result())
