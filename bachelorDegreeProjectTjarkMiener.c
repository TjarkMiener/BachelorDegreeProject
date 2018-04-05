//
//  bachelorDegreeProjectTjarkMiener.c
//
//  Created by Tjark Miener and Bruce Allen on 29.07.14.
//  Copyright (c) 2014 Tjark Miener and Bruce Allen. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define Allocate_blocksize 1024
#define Stripestheta 3001
#define Stripesphi (2*Stripestheta)
#define Gamma_max (2*M_PI/Stripesphi)
#define Empty -12345
#define End_of_first_list 1265589
#define Minimum_frequency 0.01
#define Maximum_scan 64

int num_in_stripe[Stripestheta];

//The structure skypoint have two doubles for the position and two ints for
//the cell the skypoint is in. Every skypoint has a number (int) and a frequence (double). 

struct Skypoint{
    int number;
    char id[Maximum_scan];
  
    double phi;      //0,...,2*PI
    double theta;    //0,...,PI
    
    int p_cell;
    int t_cell;

    double frequency;
    float significance;
};
struct Skypoint *list=NULL;

//function to compute the cells for each skypoint

void cellcomputing(struct Skypoint *a){
    
    if ((a->theta)<M_PI) {
        a->t_cell= (int)((a->theta)*Stripestheta/M_PI);
    } else {
        a->t_cell=Stripestheta-1;
    }
    a->p_cell=floor((a->phi)*(num_in_stripe[a->t_cell]-1)/(2.0*M_PI));
    return;
}

//calculate gamma of the skypoints

double angle(struct Skypoint *a,struct Skypoint *b) {
    double x = sin(a->theta)*sin(b->theta)*cos((a->phi)-(b->phi))+cos(a->theta)*cos(b->theta);
    if (x > 1.0) {
        x = 1.0;
    }
    return acos(x)*(180.0/M_PI);
}

//the function return if we need to check the cells
//the runnercell is always below the target cell

int close_enough(struct Skypoint *a, struct Skypoint *b){
    
    struct Skypoint targetcell_left_bottom;
    struct Skypoint targetcell_right_bottom;
    struct Skypoint runnercell_left_top;
    struct Skypoint runnercell_right_top;
   
    int p1=a->p_cell;
    int nextp1=(a->p_cell)+1;
    int nextt1=(a->t_cell)+1;
    
    int p2=b->p_cell;
    int t2=b->t_cell;
    int nextp2=(b->p_cell)+1;
    
    targetcell_right_bottom.theta=targetcell_left_bottom.theta=nextt1*M_PI/Stripestheta;
    targetcell_left_bottom.phi=2*M_PI*p1/num_in_stripe[a->t_cell];
    targetcell_right_bottom.phi=2*M_PI*nextp1/num_in_stripe[a->t_cell];
    
    runnercell_right_top.theta=runnercell_left_top.theta=t2*M_PI/Stripestheta;
    runnercell_left_top.phi=2*M_PI*p2/num_in_stripe[b->t_cell];
    runnercell_right_top.phi=2*M_PI*nextp2/num_in_stripe[b->t_cell];
    
    double gamma_1=angle(&targetcell_left_bottom,&runnercell_left_top);
    double gamma_2=angle(&targetcell_left_bottom,&runnercell_right_top);
    double gamma_3=angle(&targetcell_right_bottom,&runnercell_left_top);
    double gamma_4=angle(&targetcell_right_bottom,&runnercell_right_top);
    
    if (gamma_1<Gamma_max || gamma_2<Gamma_max || gamma_3<Gamma_max || gamma_4<Gamma_max) {
        return 1;
    } else {
        return 0;
    }
}

//dynamic allocation function

void my_Allocate(int a) {
    list=realloc(list,a*sizeof(struct Skypoint));
    if (!list) {
        printf("Error! No more RAM available, skypoints_max = %d!\n", a);
        exit(1);
    }
}

//comparefunktion for the qsort() of <stdlib.h>

int compare(const void *a, const void *b){
    const struct Skypoint *elementa = a;
    const struct Skypoint *elementb = b;
    
    if (elementa->t_cell > elementb->t_cell) return 1;
    if (elementa->t_cell < elementb->t_cell) return -1;
    
    if (elementa->p_cell > elementb->p_cell) return 1;
    if (elementa->p_cell < elementb->p_cell) return -1;
    
    return 0;
}

//print out function

void print_output(int i,int j,double g){
    
    char string1[Maximum_scan];
    char string2[Maximum_scan];
    
    strncpy(string1, list[i].id, Maximum_scan);
    strncpy(string2, list[j].id, Maximum_scan);
    
    char *ptr1;
    char *ptr2;
    
    if (list[i].significance==0.0 && list[j].significance!=0.0) {
        
        if (fabs(list[i].frequency-list[j].frequency) < Minimum_frequency){
            
            ptr1 = strtok(string1,"b");
            ptr2 = strtok(string2,"b");
            
            int result = strcmp (ptr1,ptr2);
            if (result != 0) {
                printf("%s %s %f %f %f %f %f %f %f %f\n",
                       list[j].id,list[i].id,list[j].frequency,list[i].frequency,
                       list[j].phi,list[j].theta,list[i].phi,list[i].theta,
                       list[j].significance,g);
            }
        }
        
    } else if (list[j].significance==0.0 && list[i].significance!=0.0){
        
        if (fabs(list[i].frequency-list[j].frequency) < Minimum_frequency){
            
            ptr1 = strtok(string1,"b");
            ptr2 = strtok(string2,"b");
            
            int result = strcmp (ptr1,ptr2);
            if (result != 0) {
                printf("%s %s %f %f %f %f %f %f %f %f\n",
                       list[i].id,list[j].id,list[i].frequency,list[j].frequency,
                       list[i].phi,list[i].theta,list[j].phi,list[j].theta,
                       list[i].significance,g);
            }
        }
    }
    return;
}

//main function

int main (int argc, char** argv) {
    
    int i,j,k,l,m;
    
    int skypoints_read=0;
    int skypoints_allocated=0;
    
    int target;
    int runner,runner_below_begin,runner_below_end;
    int check_point;
    int first_point_in_stripe[Stripestheta+1];
    int last_stripe;
    
    double gamma;
    
    //check if Stripestheta is odd
    if (Stripestheta%2==0) {
        printf("Error! Stripestheta = %d is not odd!\n", Stripestheta);
        exit(2);
    }
    
    //initialisation of num_in_stripe[]
    
    for (i=0; i<Stripestheta; i++) {
        num_in_stripe[i]=ceil(Stripesphi*sin(M_PI*(i+0.5)/Stripestheta));
    }
    
    //check if arguments in the command line is right
    
    if (argc!=2) {
        printf("Error! Missing command line arguments, only found %d not 1\n", argc-1);
        exit(3);
    }

    //open the .txt file
    
    FILE *input_candidate = fopen(argv[1],"r");
    
    if (!input_candidate) {
        printf("Error! File %s not found\n", argv[1]);
        exit(4);
    }
    
    while (true) {
        
        //reallocating memory if necessary
        if (skypoints_read==skypoints_allocated) {
            skypoints_allocated += Allocate_blocksize;
            my_Allocate(skypoints_allocated);
        }
        
        if (skypoints_read < End_of_first_list) {
            
            //scan the elements of the first part of the file
            if (6==(i=fscanf(input_candidate, "%d %s %lf %lf %lf %f",&list[skypoints_read].number,list[skypoints_read].id,&list[skypoints_read].phi,&list[skypoints_read].theta,&list[skypoints_read].frequency,&list[skypoints_read].significance))) {
                
                list[skypoints_read].phi=(list[skypoints_read].phi*15.0)*2.0*M_PI/360.0;
                list[skypoints_read].theta=(90.0-list[skypoints_read].theta)*2.0*M_PI/360.0;
                
                //check the data ranges, exit if error
                if (list[skypoints_read].phi < 0 || list[skypoints_read].phi > 2*M_PI || list[skypoints_read].theta < 0 || list[skypoints_read].theta > M_PI) {
                    printf("Error! Problem reading %s at line %d. Point (%f,%f) isn't in the data range!\n", argv[1],skypoints_read, list[skypoints_read].phi, list[skypoints_read].theta);
                    exit(5);
                }
                
                cellcomputing(list+skypoints_read);
                skypoints_read++;
                
            } else break;
            
        } else {
            
            //scan the elements of the secound part of the file
            if (5==(i=fscanf(input_candidate, "%d %s %lf %lf %lf",&list[skypoints_read].number,list[skypoints_read].id,&list[skypoints_read].phi,&list[skypoints_read].theta,&list[skypoints_read].frequency))) {
                
                list[skypoints_read].phi=list[skypoints_read].phi*2.0*M_PI/360.0;
                list[skypoints_read].theta=(90.0-list[skypoints_read].theta)*2.0*M_PI/360.0;
                list[skypoints_read].significance=0.0;
                
                //check the data ranges, exit if error
                if (list[skypoints_read].phi < 0 || list[skypoints_read].phi > 2*M_PI || list[skypoints_read].theta < 0 || list[skypoints_read].theta > M_PI) {
                    printf("Error! Problem reading %s at line %d. Point (%f,%f) isn't in the data range!\n", argv[1],skypoints_read, list[skypoints_read].phi, list[skypoints_read].theta);
                    exit(5);
                }
                
                cellcomputing(list+skypoints_read);
                skypoints_read++;
                
            } else break;
            
        }
        
    }
    
    if (i!=EOF) {
        fprintf(stderr, "Error! Problem reading %s at line %d (number:%d). fscanf() returned %d\n", argv[1], skypoints_read, list[skypoints_read].number,i);
        exit(6);
    }

  
    //adding shadow points
    
    k=j=0;
    for (i=0; i<skypoints_read; i++) {
        
        //reallocating memory if necessary
        if ((skypoints_read+k)==skypoints_allocated && j==1) {
            skypoints_allocated += Allocate_blocksize;
            my_Allocate(skypoints_allocated);
        }
        
        //adding shadow points on the right
        if (list[i].p_cell < 2){
         
            list[skypoints_read+k].phi=((list[i].phi)+2.0*M_PI);
            list[skypoints_read+k].theta=list[i].theta;
            list[skypoints_read+k].t_cell=list[i].t_cell;
            list[skypoints_read+k].p_cell=list[i].p_cell+num_in_stripe[list[i].t_cell];
            list[skypoints_read+k].significance=list[i].significance;
            list[skypoints_read+k].frequency=list[i].frequency;
            list[skypoints_read+k].number=list[i].number;
            strncpy(list[skypoints_read+k].id, list[i].id, Maximum_scan);
            
            k++;
            j=1;
        
        //adding shadow points on the left
        } else if (list[i].p_cell > num_in_stripe[list[i].t_cell]-3) {
         
            list[skypoints_read+k].phi=((list[i].phi)-2.0*M_PI);
            list[skypoints_read+k].theta=list[i].theta;
            list[skypoints_read+k].t_cell=list[i].t_cell;
            list[skypoints_read+k].p_cell=list[i].p_cell-num_in_stripe[list[i].t_cell];
            list[skypoints_read+k].significance=list[i].significance;
            list[skypoints_read+k].frequency=list[i].frequency;
            list[skypoints_read+k].number=list[i].number;
            strncpy(list[skypoints_read+k].id, list[i].id, Maximum_scan);
         
            k++;
            j=1;
        } else {
            j=0;
        }

    }
    
    skypoints_read += k;

    //sorting the list with the shadow points
    
    qsort (list,skypoints_read,sizeof(struct Skypoint), compare);
    
    //find the first point in the stripe
    
    for (i=0; i<=Stripestheta; i++) {
        first_point_in_stripe[i]=Empty;
    }
    
    last_stripe=first_point_in_stripe[list[1].t_cell]=0;
    for (i=1; i<skypoints_read; i++) {
        if (list[i].t_cell != last_stripe) {
            last_stripe=list[i].t_cell;
            first_point_in_stripe[list[i].t_cell]=i;
        }
    }

    //compute gamma of the points in the same cell and in the next cell
    
    for (target=0; target<skypoints_read; target++) {
      
        //only take the non shadow point
        if (list[target].phi>=0 && list[target].phi<=2.0*M_PI) {
        
            for (runner=target+1; runner<skypoints_read; runner++) {
                
                //target and runner are in the same cell
                if (list[target].p_cell==list[runner].p_cell && list[target].t_cell==list[runner].t_cell) {
                    gamma=angle(&list[target],&list[runner]);
                    print_output(target,runner,gamma);
                
                //runner are in the next cell
                } else if(list[target].p_cell==((list[runner].p_cell)-1) && list[target].t_cell==list[runner].t_cell) {
                    gamma=angle(&list[target],&list[runner]);
                    print_output(target,runner,gamma);
        
                //runner are too far away -> break
                } else break;
            }
        }
    }
    
    
    //initialization of the first point
    
    for (i=0; i<skypoints_read; i++) {
        
        //only take the non shadow point
        if (list[i].phi>=0 && list[i].phi<=2.0*M_PI) {
        
            //if the stripe below is empty, take the next point in the For-loop and check!
            if (first_point_in_stripe[((list[i].t_cell)+1)]==Empty) continue;
        
            //check if any point in the stripe below is close_enough()
            for (j=first_point_in_stripe[list[i].t_cell+1]; list[j].t_cell==(list[i].t_cell+1); j++) {
                if (close_enough(&list[i],&list[j])) {
                    runner_below_begin=j;
                    check_point=j;
                    break;
                } else check_point=Empty;
            }
            
            //if no point in the stripe below is close_enough(), take the next point in the For-loop and check
            if (check_point==Empty) continue;
        
            //check how much points in the stripe below are close_enough()
            for (k=runner_below_begin; list[k].t_cell==(list[i].t_cell+1); k++) {
                if (close_enough(&list[i],&list[k])) {
                    runner_below_end=k;
                } else break;
            }

            break;
        }
    }
    
    //if we can't initial the first point, it means that for every point there is no point in the stripe below close_enough()
    //we can skip the following part and jump with goto to end of the programm
    if (i==skypoints_read) goto end_of_program;
    
    //compute gamma of the points in the stripe below which are close_enough() to the target cell
    
    for (runner=runner_below_begin; runner<=runner_below_end; runner++) {
        gamma=angle(&list[i],&list[runner]);
        print_output(i,runner,gamma);
    }
    
    //after initialize the first point, compute the other points
    
    for (target=i+1; target<skypoints_read; target++) {
       
        //only take the non shadow point
        if (list[target].phi>=0 && list[target].phi<=2.0*M_PI) {
            
            //if the stripe below is empty, take the next point in the For-loop and check
            if (first_point_in_stripe[((list[target].t_cell)+1)]==Empty) continue;

            //the target stays in the same cell (no update of the runner_below_begin and runner_below_end)
            if (list[target].p_cell==list[target-1].p_cell && list[target].t_cell==list[target-1].t_cell) {
                
                //if no point in the stripe below is close_enough(), take the next point in the For-loop and check!
                if (check_point==Empty) continue;
                
                //compute gamma of the points in the cell below which are close_enough() to the target cell
                for (runner=runner_below_begin; runner<=runner_below_end; runner++) {
                    gamma=angle(&list[target],&list[runner]);
                    print_output(target,runner,gamma);
                }
                
            //the target moves to a different cell (update of the runner_below_begin and runner_below_end)
            } else if (list[target].p_cell!=list[target-1].p_cell && list[target].t_cell==list[target-1].t_cell) {
                
                //check if any point in the stripe below is close_enough()
                for (j=first_point_in_stripe[((list[target].t_cell)+1)]; list[j].t_cell==list[target].t_cell+1; j++) {
                    if (close_enough(&list[target],&list[j])) {
                        runner_below_begin=j;
                        check_point=j;
                        break;
                    } else check_point=Empty;
                }
                
                //if no point in the stripe below is close_enough(), take the next point in the For-loop and check!
                if (check_point==Empty) continue;
                
                //check how much points in the stripe below are close_enough()
                for (k=runner_below_begin; list[k].t_cell==(list[target].t_cell+1); k++) {
                    if (close_enough(&list[target],&list[k])) {
                        runner_below_end=k;
                    } else break;
                }
                
                //compute gamma of the points in the stripe below which are close_enough() to the target cell
                for (runner=runner_below_begin; runner<=runner_below_end; runner++) {
                    gamma=angle(&list[target],&list[runner]);
                    print_output(target,runner,gamma);
                }
            
            //the target moves to a different stripe (update of the runner_below_begin and runner_below_end)
            } else {
            
                //check if any point in the stripe below is close_enough()
                for (l=first_point_in_stripe[((list[target].t_cell)+1)]; list[l].t_cell==(list[target].t_cell+1); l++) {
                    if (close_enough(&list[target],&list[l])) {
                        runner_below_begin=l;
                        check_point=l;
                        break;
                    } else check_point=Empty;
                }
            
                //if no point in the stripe below is close_enough(), take the next point in the For-loop and check!
                if (check_point==Empty) continue;
            
                //check how much points in the stripe below are close_enough()
                for (m=runner_below_begin; list[m].t_cell==list[target].t_cell+1; m++) {
                    if (close_enough(&list[target],&list[m])) {
                        runner_below_end=m;
                    } else break;
                }
                
                //compute gamma of the points in the stripe below which are close_enough() to the target cell
                for (runner=runner_below_begin; runner<=runner_below_end; runner++) {
                    gamma=angle(&list[target],&list[runner]);
                    print_output(target,runner,gamma);
                }
            }
        }
    }
    
    end_of_program : printf("There are no more points to compare!\n");
  
    free(list);
    fclose(input_candidate);
    
    return 0;
}
