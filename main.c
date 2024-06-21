#ifndef _XOPEN_SOURCE_EXTENDED
#define _XOPEN_SOURCE_EXTENDED 1
#endif
#define PDC_WIDE

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <curses.h>

#ifdef WACS_S1
# define HAVE_WIDE 1
#else
# define HAVE_WIDE 0
#endif

#include <locale.h>

#if HAVE_WIDE
# include <wchar.h>
#endif

struct commands
{
    const char *text;
    void (*function)(WINDOW *);
};
typedef struct commands COMMAND;

#define MAXINPUTLEN 50
typedef struct inputbox{
    char input_str[MAXINPUTLEN];
    struct inputbox * above_input;
    struct inputbox * below_input;
    struct inputbox * left_input;
    struct inputbox * right_input;
    int rgroup;
    int ridx;
    int y;
    int x;
    int width;
    char type; // i : int, d : double, s : string, r : radio, q : confirm button, p : previous page, n : next page

} InputBox;

InputBox make_input(char type,int y,int x,int width,char * init_str){
    InputBox new_input = {"",NULL,NULL,NULL,NULL,0,0,y,x,width,type};
    strcpy(new_input.input_str,init_str);
    return new_input;
}

InputBox make_radio_input(int rgroup,int ridx,int y,int x,int width,char * option_str){
    InputBox new_input = {"",NULL,NULL,NULL,NULL,rgroup,ridx,y,x,width,'r'};
    strcpy(new_input.input_str,option_str);
    return new_input;
}

int input_to_int(InputBox * target){
    return atoi(target->input_str);
}

int selected_default_units[3] = {-1,-1,-1};

double input_to_double(InputBox * target,char type){
    char * unit;
    double result = strtod(target->input_str,&unit);
    switch(type)
    {
        case 'm':
            switch(selected_default_units[0])
            {
                case 0: //kg
                    break;
                case 1: //Solar Mass
                    result *= 1.98847e+30;
                    break;
                case 2: //Jupiter Mass
                    result *= 1.89813e+27;
                    break;
                case 3: //Earth Mass
                    result *= 5.9722e+24;
                    break;
                case 4: //Lunar Mass
                    result *= 7.342e+22;
                    break;

            }
            break;
        case 'p':
            switch(selected_default_units[1])
            {
                case 0: //m
                    break;
                case 1: //km
                    result *= 1000;
                    break;
                case 2: //au
                    result *= 149597870700;
                    break;
            }
            break;
        case 'v':
            switch(selected_default_units[2])
            {
                case 0: //m/s
                    break;
                case 1: //km/s
                    result *= 1000;
                    break;
            }
            break;
    }
    return result;
}

void input_to_name(InputBox * target,int idx){
    if(strcmp(target->input_str,""))
        return;
    else{
        sprintf(target->input_str,"body_%d",idx);
        return;
    }
}

void edit_input(InputBox * target){
    color_set(6,NULL);
    mvgetnstr(target->y,target->x,target->input_str,target->width);
}

void draw_general_str(InputBox * target,short c){
    color_set(c,NULL);
    mvaddnstr(target->y,target->x,target->input_str,target->width);
}

void draw_general_bg(InputBox * target,short c){
    color_set(c,NULL);
    chtype ch = ' ';
    mvaddch(target->y,target->x,ch);
    int i;
    for(i=1;i<(target->width);i++){
        addch(ch);
    }
}

void draw_input_bg(InputBox * target){
    draw_general_bg(target,4);
}

void draw_cur_input_bg(InputBox * target){
    draw_general_bg(target,5);
}

void draw_edit_input_bg(InputBox * target){
    draw_general_bg(target,6);
}

void draw_input_str(InputBox * target){
    draw_general_str(target,4);
}

void draw_cur_input_str(InputBox * target){
    draw_general_str(target,5);
}

void draw_edit_input_str(InputBox * target){
    draw_general_str(target,6);
}

void draw_input(InputBox * target){
    draw_input_bg(target);
    draw_input_str(target);
}

void draw_cur_input(InputBox * target){
    draw_cur_input_bg(target);
    draw_cur_input_str(target);
}

void draw_edit_input(InputBox * target){
    draw_edit_input_bg(target);
    draw_edit_input_str(target);
}

void h_connect_input(InputBox * left, InputBox * right){
    left->right_input = right;
    right->left_input = left;
}

void v_connect_input(InputBox * above, InputBox * below){
    above->below_input = below;
    below->above_input = above;
}

#define _USE_MATH_DEFINES
#include <math.h>

const double G = 6.67430e-11;

typedef struct vector{
	double x;
	double y;
	double z;
} Vector;

#define MAXNAMELEN 30
typedef struct celestial{
    char name[MAXNAMELEN];
    Vector pos;
	Vector vel;
	Vector acc;
	double mass;
} Celestial;

Vector makeVec(double x, double y, double z){
	Vector result;
	result.x=x;
	result.y=y;
	result.z=z;
	return result;
}

Vector scalar_multi(double sca, Vector vec){	//sca*vec
	Vector result = vec;
	result.x *= sca;
	result.y *= sca;
	result.z *= sca;
	return result;
}

Vector scalar_div(Vector vec, double sca){	//vec/sca
	return scalar_multi(1/sca,vec);
}

Vector vec_add(Vector vec1, Vector vec2){	//vec1+vec2
	Vector result = vec1;
	result.x += vec2.x;
	result.y += vec2.y;
	result.z += vec2.z;
	return result;
}

Vector vec_subtr(Vector vec1, Vector vec2){	//vec1-vec2
	Vector result = vec1;
	result.x -= vec2.x;
	result.y -= vec2.y;
	result.z -= vec2.z;
	return result;
}

double vec_mag(Vector vec){	//|vec|
	return sqrt( pow(vec.x,2) + pow(vec.y,2) + pow(vec.z,2) );
}

Vector unit_vec(Vector vec){	//unit vector of vec
	return scalar_div(vec,vec_mag(vec));
}

Vector rel_pos(Vector pos1,Vector pos2){	//position of pos2 relative to pos1
	return vec_subtr(pos2,pos1);
}

Vector acc_gravity(Celestial c1, Celestial c2){	//c1's acceleration by gravity of c2
	Vector r = rel_pos(c1.pos,c2.pos);
	return scalar_multi( (G*c2.mass)/pow(vec_mag(r),2), unit_vec(r) );
}

double kepler_law3_period(double m1, double m2, double a){
	double mu = G*(m1+m2);
	return 2*M_PI*pow(pow(a,3)/mu, 0.5);
}

double circular_orbit_speed(double m1, double m2, double r){
	double mu = G*(m1+m2);
	//printf("mu: %f\n",mu);
	return pow(mu/r,0.5);
}

void n_body_sim(WINDOW *win){
    int i,j,k,key;
    int num_celest=0;
	double max_t=0;
	double t=0;
	double delta_t=1;

	mvprintw(2,10,"[Number of Celestial bodies]");

	InputBox num_celest_input = make_input('i',3,10,6,"");
    InputBox * cur_input = &num_celest_input;
    draw_edit_input_bg(cur_input);
    echo();
    edit_input(cur_input);
    flushinp();
    noecho();
    num_celest = input_to_int(cur_input);
    erase();
    color_set(0,NULL);
    //mvprintw(2,30,"num_celest: %d",num_celest);
    //napms(2000);

    //declare default unit inputs
    InputBox mass_default_unit_inputs[5];
    InputBox pos_default_unit_inputs[3];
    InputBox vel_default_unit_inputs[2];

    //define default unit inputs
    mvprintw(2,10,"[Mass default unit options]");
    char * mass_unit_options[5] = {"kg","Solar mass (S)","Jupiter mass (J)","Earth mass (E)","Lunar mass (L)"};
    for(i=0;i<5;i++)
        mass_default_unit_inputs[i]=make_radio_input(0,i,3,10+20*i,18,mass_unit_options[i]);

    mvprintw(5,10,"[Position default unit options]");
    char * pos_unit_options[3] = {"m","km","au"};
    for(i=0;i<3;i++)
        pos_default_unit_inputs[i]=make_radio_input(1,i,6,10+6*i,4,pos_unit_options[i]);

    mvprintw(8,10,"[Velocity default unit options]");
    char * vel_unit_options[2] = {"m/s","km/s"};
    for(i=0;i<2;i++)
        vel_default_unit_inputs[i]=make_radio_input(2,i,9,10+6*i,4,vel_unit_options[i]);

    InputBox default_unit_confirm_btn = make_input('q',11,10,10,"Confirm");

    //connect all the inputs
    for(i=0;i<4;i++)
        h_connect_input(&mass_default_unit_inputs[i],&mass_default_unit_inputs[i+1]);
    for(i=0;i<2;i++)
        h_connect_input(&pos_default_unit_inputs[i],&pos_default_unit_inputs[i+1]);
    h_connect_input(&vel_default_unit_inputs[0],&vel_default_unit_inputs[1]);

    for(i=0;i<3;i++)
        v_connect_input(&mass_default_unit_inputs[i],&pos_default_unit_inputs[i]);
    for(i=3;i<5;i++)
        mass_default_unit_inputs[i].below_input = &pos_default_unit_inputs[2];
    for(i=0;i<2;i++)
        v_connect_input(&pos_default_unit_inputs[i],&vel_default_unit_inputs[i]);
    pos_default_unit_inputs[2].below_input = &vel_default_unit_inputs[1];
    for(i=1;i>=0;i--)
        v_connect_input(&vel_default_unit_inputs[i],&default_unit_confirm_btn);

    //start default unit selection page
    for(i=0;i<3;i++)
        selected_default_units[i] = -1; // <- this array is a global variable defined before input_to_double function
    cur_input = &mass_default_unit_inputs[0];
    bool repeat = 1;
    while(repeat)
    {
        for(i=0;i<5;i++)
            draw_input(&mass_default_unit_inputs[i]);

        for(i=0;i<3;i++)
            draw_input(&pos_default_unit_inputs[i]);

        for(i=0;i<2;i++)
            draw_input(&vel_default_unit_inputs[i]);

        draw_input(&default_unit_confirm_btn);

        draw_cur_input(cur_input);

        if(selected_default_units[0]>=0)
            draw_edit_input(&mass_default_unit_inputs[selected_default_units[0]]);
        if(selected_default_units[1]>=0)
            draw_edit_input(&pos_default_unit_inputs[selected_default_units[1]]);
        if(selected_default_units[2]>=0)
            draw_edit_input(&vel_default_unit_inputs[selected_default_units[2]]);

        noecho();
        keypad(stdscr, TRUE);
        raw();
        key = getch();
        switch(key)
        {
            case KEY_UP:
                if(cur_input->above_input!=NULL)
                    cur_input = cur_input->above_input;
                break;
            case KEY_DOWN:
                if(cur_input->below_input!=NULL)
                    cur_input = cur_input->below_input;
                break;
            case KEY_LEFT:
                if(cur_input->left_input!=NULL)
                    cur_input = cur_input->left_input;
                break;
            case KEY_RIGHT:
                if(cur_input->right_input!=NULL)
                    cur_input = cur_input->right_input;
                break;

            case 10:
            case 13:
            case KEY_ENTER:
                switch(cur_input->type)
                {
                    case 'r':
                        selected_default_units[cur_input->rgroup] = cur_input->ridx;
                        break;
                    case 'q':
                        repeat=0;
                }
            default:
                break;
        }
    }

    erase();
    color_set(0,NULL);
    //mvprintw(3,30,"%d // %d / %d / %d",num_celest,selected_default_units[0],selected_default_units[1],selected_default_units[2]);
    //refresh();
    //napms(5000);


    //declare inputs
    InputBox * name_inputs = (InputBox *)malloc(sizeof(InputBox)*num_celest);    //InputBox name_inputs[num_celest];
    InputBox * mass_inputs = (InputBox *)malloc(sizeof(InputBox)*num_celest);    //InputBox mass_inputs[num_celest];

    //InputBox init_pos_inputs[num_celest][3];
    //InputBox init_vel_inputs[num_celest][3];
    InputBox ** init_pos_inputs = (InputBox **)malloc(sizeof(InputBox *) * num_celest);
    InputBox ** init_vel_inputs = (InputBox **)malloc(sizeof(InputBox *) * num_celest);
    for(i=0;i<num_celest;i++){
        init_pos_inputs[i] = (InputBox *)malloc(sizeof(InputBox)*3);
        init_vel_inputs[i] = (InputBox *)malloc(sizeof(InputBox)*3);
    }
    InputBox * next_buttons = (InputBox *)malloc(sizeof(InputBox)*num_celest);     //InputBox next_buttons[num_celest];
    InputBox * prev_buttons = (InputBox *)malloc(sizeof(InputBox)*num_celest);     //InputBox prev_buttons[num_celest];

    //define inputs
    #define START_LINE 4
    #define LINE_PAD 3
    for(i=0;i<num_celest;i++){
        int input_y = START_LINE;
        if(!i)
            mvaddstr(input_y-1,10,"[Name]");
        name_inputs[i] = make_input('s',input_y,10,30,"");
        input_y += LINE_PAD;
        if(!i)
            mvaddstr(input_y-1,10,"[Mass]");
        mass_inputs[i] = make_input('d',input_y,10,33,"");
        input_y += LINE_PAD;
        if(!i){
            mvaddstr(input_y-1,10,"[Initial Position] (x,y,z)");
            mvaddstr(input_y-1+LINE_PAD,10,"[Initial Velocity] (x,y,z)");
        }
        for(j=0;j<3;j++){
            init_pos_inputs[i][j] = make_input('d',input_y,10+22*j,20,"");
            init_vel_inputs[i][j] = make_input('d',input_y+LINE_PAD,10+22*j,20,"");
        }
        next_buttons[i]=make_input('n',START_LINE+2*LINE_PAD,76,2,"->");
        prev_buttons[i]=make_input('p',START_LINE+2*LINE_PAD,6,2,"<-");
    }

    //connect inputs
    for(i=0;i<num_celest;i++){
        h_connect_input(prev_buttons+i,&init_pos_inputs[i][0]);
        h_connect_input(&init_pos_inputs[i][2],next_buttons+i);

        name_inputs[i].left_input = prev_buttons+i;
        mass_inputs[i].left_input = prev_buttons+i;
        name_inputs[i].right_input = next_buttons+i;
        mass_inputs[i].right_input = next_buttons+i;
        init_vel_inputs[i][0].left_input = prev_buttons+i;
        init_vel_inputs[i][2].right_input = next_buttons+i;

        for(j=0;j<2;j++){
            h_connect_input(&init_pos_inputs[i][j],&init_pos_inputs[i][j+1]);
            h_connect_input(&init_vel_inputs[i][j],&init_vel_inputs[i][j+1]);
        }

        v_connect_input(name_inputs+i,mass_inputs+i);
        v_connect_input(mass_inputs+i,&init_pos_inputs[i][0]);
        for(j=1;j<3;j++){
            init_pos_inputs[i][j].above_input=mass_inputs+i;
        }
        for(j=0;j<3;j++){
            v_connect_input(&init_pos_inputs[i][j],&init_vel_inputs[i][j]);
        }
    }

    //start the page
    i=0;
    cur_input = &name_inputs[0];
    while(i<num_celest)
    {
        draw_input(name_inputs+i);
        draw_input(mass_inputs+i);
        for(j=0;j<3;j++){
            draw_input(&init_pos_inputs[i][j]);
            draw_input(&init_vel_inputs[i][j]);
        }
        draw_input(prev_buttons+i);
        draw_input(next_buttons+i);
        draw_cur_input(cur_input);

        noecho();
        keypad(stdscr, TRUE);
        raw();
        key = getch();
        switch(key)
        {
            case KEY_UP:
                if(cur_input->above_input!=NULL)
                    cur_input = cur_input->above_input;
                break;
            case KEY_DOWN:
                if(cur_input->below_input!=NULL)
                    cur_input = cur_input->below_input;
                break;
            case KEY_LEFT:
                if(cur_input->left_input!=NULL)
                    cur_input = cur_input->left_input;
                break;
            case KEY_RIGHT:
                if(cur_input->right_input!=NULL)
                    cur_input = cur_input->right_input;
                break;

            case 10:
            case 13:
            case KEY_ENTER:
                switch(cur_input->type)
                {
                    case 's':
                    case 'd':
                        draw_edit_input_bg(cur_input);
                        echo();
                        edit_input(cur_input);
                        flushinp();
                        noecho();
                        break;
                    case 'p':
                        if(i){
                            i--;
                            cur_input = name_inputs + i;
                        }
                        else
                            return;
                        break;
                    case 'n':
                        i++;
                        if(i<num_celest)
                            cur_input = name_inputs + i;
                        break;
                }
            default:
                break;
        }
    }

    Celestial * celest = (Celestial *)malloc(sizeof(Celestial)*num_celest);

    for(i=0;i<num_celest;i++){
        input_to_name(name_inputs+i,i);
        strcpy(celest[i].name,name_inputs[i].input_str);

        celest[i].pos = makeVec(input_to_double(&init_pos_inputs[i][0],'p'),
                                input_to_double(&init_pos_inputs[i][1],'p'),
                                input_to_double(&init_pos_inputs[i][2],'p'));

        celest[i].vel = makeVec(input_to_double(&init_vel_inputs[i][0],'v'),
                                input_to_double(&init_vel_inputs[i][1],'v'),
                                input_to_double(&init_vel_inputs[i][2],'v'));

        celest[i].acc = makeVec(0,0,0);
        celest[i].mass = input_to_double(mass_inputs+i,'m');
    }

    free(name_inputs);  free(mass_inputs);
    free(next_buttons); free(prev_buttons);
    for(i=0;i<num_celest;i++){
        free(init_pos_inputs[i]);
        free(init_vel_inputs[i]);
    }
    free(init_pos_inputs);
    free(init_vel_inputs);


    /*
    InputBox max_t_input;
    InputBox delta_t_input;
    InputBox file_name_input;
    */

}

	/*
	scanw("%d",&num_celest);

	Celestial celest[num_celest];
	for(i=0;i<num_celest;i++){

		mvprintw(3+i,30,"celest[%d]:\n",i);

		mvprintw(3+i+1,30,"\tmass: ");
		refresh();
		scanw("%lf", &celest[i].mass);

		mvprintw(3+i+2,30,"\tposition (x,y,z): ");
		refresh();
		scanw("%lf %lf %lf", &celest[i].pos.x, &celest[i].pos.y, &celest[i].pos.z);

		mvprintw(3+i+3,30,"\tvelocity (x,y,z): ");
		refresh();
		scanw("%lf %lf %lf", &celest[i].vel.x, &celest[i].vel.y, &celest[i].vel.z);

	}

	erase();

	for(i=0;i<num_celest;i++){

		mvprintw(i+1,30,"celest[%d]:\n",i);

		mvprintw(i+2,30,"\tmass: ");
		mvprintw(i+2,30,"%lf\n", celest[i].mass);

		mvprintw(i+3,30,"\tposition (x,y,z): ");
		mvprintw(i+3,30,"%lf %lf %lf\n", celest[i].pos.x, celest[i].pos.y, celest[i].pos.z);

		mvprintw(i+4,30,"\tvelocity (x,y,z): ");
		mvprintw(i+4,30,"%lf %lf %lf\n", celest[i].vel.x, celest[i].vel.y, celest[i].vel.z);

	}

	getch();
	erase();

	mvprintw(1,30,"max_t: ");
	refresh();
	scanw("%lf",&max_t);

	mvprintw(2,30,"delta_t: ");
	refresh();
	scanw("%lf",&delta_t);

	char file_name[50]={0};
	mvprintw(3,30,"file_name: ");
	refresh();
	scanw("%s",file_name);

	FILE * fp = fopen(file_name,"wt");

	for(i=0;i<num_celest;i++){
		fprintf(fp,"%f\t",celest[i].mass);
	}
	fprintf(fp,"\n");

	for(i=0;i<num_celest;i++){
		fprintf(fp,"%f\t%f\t%f\n", celest[i].vel.x, celest[i].vel.y, celest[i].vel.z);
	}
	fprintf(fp,"\n");

	double half_delta_t = delta_t/2;

	for(t=0; t<=max_t; t+=delta_t){
		for(i=0; i<num_celest; i++){

			celest[i].acc.x=0;
			celest[i].acc.y=0;
			celest[i].acc.z=0;

			for(j=0;j<num_celest;j++){
				if(j==i)
					continue;
				celest[i].acc = vec_add( celest[i].acc, acc_gravity(celest[i],celest[j]) );	//calculate the sum of forces
			}

			fprintf(fp,"%f\t%d\t%f\t%f\t%f\n",t,i,celest[i].pos.x,celest[i].pos.y,celest[i].pos.z);	//print result

			celest[i].vel = vec_add( celest[i].vel, scalar_multi( half_delta_t, celest[i].acc ) );	//v=v0+adt

		}

		for(k=0; k<num_celest; k++){
			celest[k].pos = vec_add( celest[k].pos, scalar_multi( delta_t, celest[k].vel ) );	//x=x0+vdt	| x=x0+v0dt+(1/2)a(dt)^2 ??
		}


		for(i=0; i<num_celest; i++){

			celest[i].acc.x=0;
			celest[i].acc.y=0;
			celest[i].acc.z=0;

			for(j=0;j<num_celest;j++){
				if(j==i)
					continue;
				celest[i].acc = vec_add( celest[i].acc, acc_gravity(celest[i],celest[j]) );	//calculate the sum of forces
			}

			celest[i].vel = vec_add( celest[i].vel, scalar_multi( half_delta_t, celest[i].acc ) );	//v=v0+adt

		}

	}
}
*/
void cr3bp_sim(WINDOW *win){

}

#define INTENTIONALLY_UNUSED_PARAMETER( param) (void)(param)
#define MAX_OPTIONS (2)

COMMAND command[MAX_OPTIONS] =
{
    {"N-body", n_body_sim},
    {"CR3BP", cr3bp_sim},
};


int width, height;
static bool report_mouse_movement = FALSE;
static SCREEN *screen_pointer;

int initTest(WINDOW **win, int argc, char *argv[])
{
    INTENTIONALLY_UNUSED_PARAMETER( argv);
    INTENTIONALLY_UNUSED_PARAMETER( argc);
    screen_pointer = newterm(NULL, stdout, stdin);

    if (has_colors())
        start_color();

    /* Create a drawing window */

    width  = 60;
    height = 13;

    *win = newwin(height, width, (LINES - height) / 2, (COLS - width) / 2);

    if (*win == NULL)
    {
        endwin();
        return 1;
    }

    return 0;
}

#define UNICODE_BBLOCK      0x2584
#define UNICODE_UBLOCK      0x2580
#define UNICODE_FBLOCK      0x2588

void loading_screen(FILE *f)
{

    int y0 = 1;
    int x0 = 30;

    short c;
    for(c=1;c<=3;c++){
        int x,y;
        fread((void *)(&y),sizeof(int),1,f);
        fread((void *)(&x),sizeof(int),1,f);

        int y_flip = 51-y+1;
        int y_line = ((y_flip-1)/2)+1;
        bool ub = y_flip%2;  // 0 when BBLOCK, 1 when UBLOCK

        chtype ch,prev_ch;
        prev_ch = mvinch(y0+y_line,x0+x);

        if(ub){
            if(prev_ch==(UNICODE_BBLOCK | COLOR_PAIR(c)) || prev_ch==(UNICODE_FBLOCK | COLOR_PAIR(c)))
                ch = UNICODE_FBLOCK | COLOR_PAIR(c);
            else
                ch = UNICODE_UBLOCK | COLOR_PAIR(c);
        }
        else{
            if(prev_ch==(UNICODE_UBLOCK | COLOR_PAIR(c)) || prev_ch==(UNICODE_FBLOCK | COLOR_PAIR(c)))
                ch = UNICODE_FBLOCK | COLOR_PAIR(c);
            else
                ch = UNICODE_BBLOCK | COLOR_PAIR(c);
        }

        mvaddch(y0+y_line,x0+x,ch);
    }

    refresh();

}

void display_menu(int old_option, int new_option)
{
    int lmarg = (COLS - 14) / 2,
        tmarg = (LINES - (MAX_OPTIONS + 2)) / 2;

    if (old_option == -1)
    {
        int i;

        attrset(A_BOLD | A_REVERSE);
        mvprintw(tmarg-2, lmarg-10,"N-body & CR3BP Simulation Program");
        attrset(A_NORMAL);

        for (i = 0; i < MAX_OPTIONS; i++)
            mvaddstr(tmarg + i, lmarg, command[i].text);
    }
    else
        mvaddstr(tmarg + old_option, lmarg, command[old_option].text);

    attrset(A_REVERSE);
    mvaddstr(tmarg + new_option, lmarg, command[new_option].text);
    attrset(A_NORMAL);

    mvaddstr(tmarg + MAX_OPTIONS + 2, lmarg - 23,
             "Use Up and Down Arrows to select - Enter to run - Q to quit");
    refresh();
}

int main(int argc, char *argv[])
{
    WINDOW *win;
    int key, old_option = -1, new_option = 0, i;
    bool quit = FALSE;

#ifdef _WIN32
    setlocale(LC_ALL, ".utf8");
#endif

#ifdef __PDCURSESMOD__
#ifdef PDC_VER_MAJOR   /* so far only seen in 4.0+ */
    PDC_set_resize_limits( 20, 50, 70, 200);
#endif
#endif

    if (initTest(&win, argc, argv))
        return 1;

    for( i = 1; i < argc; i++)
        if( argv[i][0] == '-')
            switch( argv[i][1])
            {
                case 'l': case 'L':
                    setlocale( LC_ALL, argv[i] + 2);
                    break;
#ifdef __PDCURSESMOD__
                case 'b': case 'B':
                    PDC_set_blink( TRUE);
                    break;
                case 'm': case 'M':
                    PDC_return_key_modifiers( TRUE);
                    break;
                case 't':
                    traceon( );
                    break;
#ifdef PDC_VER_MAJOR   /* so far only seen in 4.0+ */
                case 'r':     /* allow user-resizable windows */
                    {
                        int min_lines, max_lines, min_cols, max_cols;

                        if( sscanf( argv[i] + 2, "%d,%d,%d,%d",
                                       &min_lines, &max_lines,
                                       &min_cols, &max_cols) == 4)
                            PDC_set_resize_limits( min_lines, max_lines,
                                                   min_cols, max_cols);
                    }
                    break;
#endif
#endif
                case 'z':
                    report_mouse_movement = TRUE;
                    break;
                default:
                    break;
            }

    erase();

    PDC_set_blink(TRUE);
    attrset(A_BLINK);
    mvaddstr(24, 52, "Press ENTER to start");

    FILE *file = fopen("planar_3body_tui.dat", "rb");
    int num_of_bodies = 0;
    int time_len = 0;
    fread((void *)(&num_of_bodies), sizeof(int), 1, file);
    fread((void *)(&time_len),sizeof(int),1,file);
    attrset(A_BOLD | A_REVERSE);
    mvprintw(1, 45,"N-body & CR3BP Simulation Program");
    attrset(A_NORMAL);
    //mvprintw(31, 30,"num_of_boies: %d, time_len: %d",num_of_bodies,time_len);

    //start_color was already called by initTest function
    init_pair(0,COLOR_WHITE,COLOR_BLACK);
    init_pair(1,COLOR_BLUE,COLOR_BLACK);
    init_pair(2,COLOR_GREEN,COLOR_BLACK);
    init_pair(3,COLOR_YELLOW,COLOR_BLACK);
    init_pair(4,COLOR_BLACK,COLOR_WHITE);
    init_pair(5,COLOR_WHITE,COLOR_YELLOW);
    init_pair(6,COLOR_WHITE,COLOR_BLUE);

    i=0;
    while(1)
    {
        loading_screen(file);
        i++;
        if(i>=352){
            napms(3000);
            break;
        }
        timeout(0);
        if (getch()==10)
            break;
        napms(15);
    }
    erase();
    display_menu(old_option, new_option);

    while (1)
    {
        noecho();
        keypad(stdscr, TRUE);
        raw();

        key = getch();

        switch(key)
        {
        case 10:
        case 13:
        case KEY_ENTER:
            old_option = -1;
            erase();
            refresh();
            (*command[new_option].function)(win);
            erase();
            display_menu(old_option, new_option);
            break;

        case KEY_PPAGE:
        case KEY_HOME:
            old_option = new_option;
            new_option = 0;
            display_menu(old_option, new_option);
            break;

        case KEY_NPAGE:
        case KEY_END:
            old_option = new_option;
            new_option = MAX_OPTIONS - 1;
            display_menu(old_option, new_option);
            break;

        case KEY_UP:
            old_option = new_option;
            new_option = (new_option == 0) ?
                new_option : new_option - 1;
            display_menu(old_option, new_option);
            break;

        case KEY_DOWN:
            old_option = new_option;
            new_option = (new_option == MAX_OPTIONS - 1) ?
                new_option : new_option + 1;
            display_menu(old_option, new_option);
            break;
#ifdef KEY_RESIZE
        case KEY_RESIZE:
# ifdef PDCURSES
            resize_term(0, 0);
# endif
            old_option = -1;
            erase();
            display_menu(old_option, new_option);
            break;
#endif
        case 'Q':
        case 'q':
            quit = TRUE;
        }

        if (quit == TRUE)
            break;
    }

    delwin(win);
    endwin();
                            /* Not really needed,  but ensures Valgrind  */
    delscreen( screen_pointer);          /* says all memory was freed */
    return 0;
}

