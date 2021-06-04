# 2021 MBL Physiology Bioinformatics Bootcamp

Welcome to the bioinformatics bootcamp!

For this course we will be using the XSEDE platform [JetStream](https://use.jetstream-cloud.org/application/projects). JetStream allows us all to be using the same "computer" and gives us more storage space and memory to carry out some bioinformatic commands.

## Getting started on JetStream

Let's get started by logging into [JetStream](https://use.jetstream-cloud.org/application/projects). Once you are logged in, you will want to go to "Projects" and select "Create New Project":

![Alt text](bash-genomics/images/js-login.png)

Create a new project and give it a name. Then, click on your project and click on "New" --> "Instance":

![Alt text](bash-genomics/images/js-new-inst.png)

We will be using an "instance" or computer set up that has been built specifically for genomics. It already has a lot of the programs we will be using pre installed. To select it you should search for "genomics":

![Alt text](bash-genomics/images/js-genomics.png)

We will be using an `m1.medium` machine that has 6 CPUs, 16Gb of RAM, and 60Gb of storage. Make sure you select that from the instance size drop down menu. There is no need to change any of the other values:


![Alt text](bash-genomics/images/js-launch.png)

The machines take a bit to start up so while they start up we are going to do some introductions and talk a bit more broadly about bioinformatics.

Once your allocation starts up we will be working within the Web Shell. To open up the terminal click on `>_Open Web Shell`:


![Alt text](bash-genomics/images/js-openws.png)

Note, if you are already comfortable with command line and want to work locally within your own terminal you may do so by setting up an ssh key. You can ask an instructor how to do that after today's class.

This should open the Web Shell in a new tab of your browser. You should see something like this:


![Alt text](bash-genomics/images/atmosphere.png)

Now we are ready to get going! Welcome to the shell.

The first thing we are going to do is download the materials we will be using from `GitHub`. To do this, copy and paste the command below into your terminal (more on how to do that in person):

```
git clone https://github.com/MBL-Physiology-Bioinformatics/2021-Bioinformatics-Tutorial-Materials.git
```

We will also run two commands to set up conda:

```
sudo su -
source /etc/bashrc
```

Now we are ready to get learning BASH! We will be working off of the [Introduction to Command Line for Genomics](https://datacarpentry.org/shell-genomics/) from Data Carpentry. Feel free to follow along there as well!

## Day 2 set up!

Today, we are going to be installing a few programs with `anaconda` ([`conda`](https://astrobiomike.github.io/unix/conda-intro)). First off, we need to get things set up within our terminal. 

First, we will run some things specific to this particular web instance: 

```
sudo su -
source /etc/bashrc
```

Then, we are going to iniate our `conda` to work with `bash`:

```
conda init bash
```

Now, we are going to re-start the shell (just refresh the broswer window) or re start your terminal connection. 

We are now ready to get going! Navigate to the `2021-Bioinformatics-Tutorial-Materials/` folder where we were working yesterday. We have made some updates to the materials today (including what you are reading now). We are going to use `git` to update our files. `git` is terribly useful for things like this. If you want to know more about how to use git for science and collaboration I encourage you to check out [this tutorial](https://github.com/swcarpentry). 

To update the tutorial materials type:

```
git pull 
```

Now, finally, we are going to run an install script in "the background". Let's take a look at what is in `setup-shell.sh`. 

```
bash setup-shell.sh > conda.log &
```

And! Onward and upward! This will run in the background while we do other things in our terminal :) 

## End of the day

We will be using the same instance each day for the course. We have a joint allocation through the XSEDE that gets charged per unit of compute per time. So, we want to make sure that we conserve our compute allocation for active use. At the end of the day it is important that everyone suspend their instance. You can stop your instance at the end of the day.  



![Alt text](bash-genomics/images/js-suspend.png)
