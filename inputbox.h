#define MAXINPUTLEN 51
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

InputBox make_input(char type,int y,int x,int width,char * init_str);

InputBox make_radio_input(int rgroup,int ridx,int y,int x,int width,char * option_str);

int input_to_int(InputBox * target);

double input_to_double(InputBox * target,char type);

void input_to_name(InputBox * target,int idx);

void edit_input(InputBox * target);

void draw_general_str(InputBox * target,short c);

void draw_general_bg(InputBox * target,short c);

void draw_input_bg(InputBox * target);

void draw_cur_input_bg(InputBox * target);

void draw_edit_input_bg(InputBox * target);

void draw_input_str(InputBox * target);

void draw_cur_input_str(InputBox * target);

void draw_edit_input_str(InputBox * target);

void draw_input(InputBox * target);

void draw_cur_input(InputBox * target);

void draw_edit_input(InputBox * target);

void h_connect_input(InputBox * left, InputBox * right);

void v_connect_input(InputBox * above, InputBox * below);
