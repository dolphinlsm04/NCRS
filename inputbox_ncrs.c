#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <locale.h>
#include <ncurses.h>
#include "inputbox.h"

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

char * mass_unit_options[5] = {"kg","Solar mass (MS)","Jupiter mass (MJ)","Earth mass (ME)","Lunar mass (ML)"};
char * pos_unit_options[3] = {"m","km","au"};
char * vel_unit_options[2] = {"m/s","km/s"};

signed char selected_default_units[3] = {-1,-1,-1};

double input_to_double(InputBox * target,char type){
    char * unit_p;
    double result = strtod(target->input_str,&unit_p);

    if(unit_p[0]>='A' && unit_p[0]<='z'){   //handle non-default units

        char unit[10] = {0};
        strncpy(unit,unit_p,9);
        int i;
        for(i=0;i<10;i++){  //make unit string lowercase

            if(unit[i]=='\0')
                break;

            unit[i] = tolower(unit[i]);
        }

        switch(type){
            case 'M':
                if(!strcmp(unit,"kg")){
                    return result;
                }
                else if(!strcmp(unit,"ms")){
                    result *= 1.98847e+30;
                    return result;
                }
                else if(!strcmp(unit,"mj")){
                    result *= 1.89813e+27;
                    return result;
                }
                else if(!strcmp(unit,"me")){
                    result *= 5.9722e+24;
                    return result;
                }
                else if(!strcmp(unit,"ml")){
                    result *= 7.342e+22;
                    return result;
                }
                else{
                    return result;
                }
                break;
            case 'p':
                if(!strcmp(unit,"m")){
                    return result;
                }
                else if(!strcmp(unit,"km")){
                    result *= 1000;
                    return result;
                }
                else if(!strcmp(unit,"au")){
                    result *= 149597870700;
                    return result;
                }
                else{
                    return result;
                }
                break;
            case 'v':
                if(!strcmp(unit,"m/s")){
                    return result;
                }
                else if(!strcmp(unit,"km/s")){
                    result *= 1000;
                    return result;
                }
                else{
                    return result;
                }
                break;

        }
    }

    switch(type)
    {
        case 'M':
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
        case 's':
            break;
        case 'm':
            result *= 60;
            break;
        case 'h':
            result *= 3600;
            break;
        case 'd':
            result *= 86400;
            break;
        case 'y':
            result *= 31557600; //Julian year - IAU standard
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
    flushinp();
    noraw();
    echo();
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
