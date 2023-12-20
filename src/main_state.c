//3tiles
//primer koji koristi razne tajlove
//tj razne slicice velicine 64x64 piksela (videt u folderu res/tiles)
//ovo mozemo da napravimo preko setova tajlova
//kod nas u primeru je namesteno da preko plusa mozemo da menjamo selektovani tajl
//
//mi  svaki od ovih tajlova imamo sacuvano kao odvojenu sliku
//dok u nekim malo efikassnijim sistemima ovo nece biti odvojene slike
//cesto cemo imati tajl setove u okviru jedne slike, jednostavno poredjane, i onda cemo da koristimo samo odredjeni deo te slike
//
//kod nas u primeru je ogranicen broj tajlova, a svaki tajl je 64x64 piksela, osim drveca, koje je nesto vece po visini 64x92
//mi cemo reci da radimo sa klasicnom velicinom, pa sve sto je vece ce biti poravnato kao da je 64x64 i gornji deo ce da viri iznad
//e sad zasto smo stavili onu cudnu roze boju, to je zato sto kada budemo iscrtavali tajlove
//kad god naidjemo na tu specificnu roze boju, nju cemo samo da ignorisemo, a sve ostalo cemo da kopiramo
//a roza se ne preslikava, nego ostaje ono sto je bilo u pozadini
//ovo je vrlo vrlo jednostavan primer transparentnosti,tj providnosti
//tj samo koristimo neku key boju, tj neku odredjenu boju, bas specificna
//u rafgl je stavljeno da je key boja (255,0,254)
//
//imamo i odredjenu matricu koja ce da predstavlja, nas svet, a svet je sastavljen iz tajlova
//tj imacemo matricu intova gde ce da stoji za koju poziciju, npr 0,0, se nalazi neki tajl, npr tajl 4
//sa plus menjamo vrednost tajla koji se prikazuje
//
//mozemo da pomeramo kameru, tj ne prikazujemo sve moguce tajlove odmah, nego pomeranjem kamere
//dobijemo prikaz i onoga sto je van ekrana, ali se to ne iscrtava dok se kamera ne pomeri na to
//
//da bismo lakse mogli da odredimo na kojoj poziciji u matrici tajlova se nalazi mis
//bitno nam je da poravno mapu tajlova sa delom koji se prikazuje od 0,0 da nam bude u gornjem levom cosku
//onda samo mozemo da podelimo poziciju misa sa sirinom i visinom tajla i da dobijemo u kom polju se nalazimo
//ali probleem nam nastaje kada prikaz i mapa nisu poravnati,zato nam je bitno gde se nalazi kamera
//a to nam je definisano preko ugla
//nama je bitno gde se nalazi 0,0 pozicija mape u odnosu na prikazan deo, tj za koliko x je pomerena mapa u levo
//i za koliko y je spusteno na dole
//i onda samo na te negativne pomeraje dodamo poziciju misa, i onda znamo u kojim pikselima u okviru mape se nalazi mis,
//tako da onda mozemo da saznamo i poziciju tog polja u mapi na kojem stoji mis
//ako je pomerano sa desne i donje strane, onda cemo da saberemo te pomeraje analogno
//tj samo da se ponasamo kao  da dovlacimo taj cosak mape u cosak odakle se prikazuje portview
//tako da ustvari racunamo (kamera_pos_x + mouse_pos_x)/vel_tajla
//kod nas u kodu je to selected_x i selected_y koje se nalaze u okviru main_state_updatea()
//i ssamo cemo imati te vrednosti jos i klempovane jer ne zelimo da izadjemo iz world sizea
#include <main_state.h>
#include <glad/glad.h>
#include <math.h>
#define GRAVITY 100.0
#define JUMP_VELOCITY 200.0


#include <rafgl.h>

#include <game_constants.h>

static rafgl_raster_t sea;
static rafgl_raster_t upscaled_sea;
static rafgl_raster_t ripple_prev, ripple_tmp;
static rafgl_raster_t raster, raster2;
static rafgl_raster_t checker;
static rafgl_raster_t cave_instance1, cave_instance1_scaled;
static rafgl_raster_t cave_instance2, cave_instance2_scaled;
static rafgl_raster_t tmpFisheye;

static rafgl_texture_t texture;

static rafgl_spritesheet_t hero;

static rafgl_spritesheet_t hero1;

//17 razlicitih tajlova i 17 rastera koji ce da predstavljaju te tajlove i koje cemo da ucitamo
#define NUMBER_OF_TILES 20
rafgl_raster_t tiles[NUMBER_OF_TILES];

//velicina matrice, koliki da nam bude svet i ta matrica
#define WORLD_SIZE 128
int tile_world[WORLD_SIZE][WORLD_SIZE];

#define TILE_SIZE 64

#define MAX_LIGHT_SOURCES

static int raster_width = RASTER_WIDTH, raster_height = RASTER_HEIGHT;

int ripple_initialized;
//pozicija karaktera na pocetku
int hero_pos_x = RASTER_WIDTH / 2 + 30;
int hero_pos_y = RASTER_HEIGHT / 2 + 30;

int hero_pos_cave_x = RASTER_WIDTH / 8 + 10;
int hero_pos_cave_y = RASTER_HEIGHT - 3 * RASTER_HEIGHT/ 8;
int ground_level;

int fisheye_cave_x;
int fisheye_cave_y;
int number_of_frames;
int fisheye_entrance;
int fisheye_number_of_frames;
int turned_on_lights;
int hero_on_ground;
int last_clicked;
int hero_vertical_velocity;

int is_cave;

//za skrinsotovanje
static char save_file[256];
int save_file_no = 0;

int time;
int frames;
int x_pos_mouse;
int y_pos_mouse;

//kamera x i kamera y je postavljeno na 600,600 na pocetku
int camx = 0, camy = 0;

int selected_x, selected_y;

int hero_is_dead = 0;

/*rafgl_pixel_rgb_t **curr;
rafgl_pixel_rgb_t **prev;
rafgl_pixel_rgb_t **temp;

void init_matrices_ripple() {
    int width = upscaled_sea.width;
    int height = upscaled_sea.height;

    temp = (int **)malloc(height * sizeof(int *));
    curr = (int **)malloc(height * sizeof(int *));
    prev = (int **)malloc(height * sizeof(int *));

    if (!temp || !curr || !prev) {
        return;
    }

    for (int i = 0; i < height; i++) {
        temp[i] = (int *)calloc(width, sizeof(int)); // calloc initializes the memory to zero
        curr[i] = (int *)calloc(width, sizeof(int));  // calloc initializes the memory to zero
        prev[i] = (int *)calloc(width, sizeof(int));  // calloc initializes the memory to zero

        if (!temp[i] || !curr[i] || !prev[i]) {
            return;
        }
    }
}

void water_ripple_effect(int mouse_pos_x, int mouse_pos_y) {
    int width = upscaled_sea.width;
    int height = upscaled_sea.height;

    /*if (mouse_pos_x < 0 || mouse_pos_x >= width || mouse_pos_y < 0 || mouse_pos_y >= height) {
        return; // Mouse click is outside the raster, so don't do anything.
    }


    const int ripple_radius = 3; // Radius of the effect.
    const int disturbance_value = 512; // Intensity of the disturbance.

    for (int dy = -ripple_radius; dy <= ripple_radius; dy++) {
        for (int dx = -ripple_radius; dx <= ripple_radius; dx++) {
            int x = mouse_pos_x + dx;
            int y = mouse_pos_y + dy;
            // Check if within bounds before applying the disturbance.
            if (dx * dx + dy * dy <= ripple_radius * ripple_radius &&
                x >= 0 && x < width && y >= 0 && y < height) {
                prev[y][x] = disturbance_value; // Access as [y][x] since it's row-major
            }
        }
    }

    float damping = 0.9;
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            // Compute the new value for the current pixel.
            int sumr = prev[y - 1][x].r + prev[y + 1][x].r + prev[y][x - 1].r + prev[y][x + 1].r;
            int sumg = prev[y - 1][x].g + prev[y + 1][x].g + prev[y][x - 1].g + prev[y][x + 1].g;
            int sumb = prev[y - 1][x].b + prev[y + 1][x].b + prev[y][x - 1].b + prev[y][x + 1].b;
            curr[y][x] = (sum / 2) - curr[y][x];
            curr[y][x].r = (int)(curr[y][x].r * damping);
            curr[y][x].g = (int)(curr[y][x].g * damping);
            curr[y][x].b = (int)(curr[y][x].b * damping);


            int value = curr[y][x];
            value = value > 255 ? 255 : value < 0 ? 0 : value;
            rafgl_pixel_rgb_t new_pixel = { .rgba = (value << 16) | (value << 8) | value }; // Assuming the format is 0xRRGGBB
            upscaled_sea.data[y * width + x] = new_pixel; // Access as [y * width + x] for row-major order
        }
    }

    int **swap = curr;
    curr = prev;
    prev = swap;
}*/

uint32_t safe_rgba_operation(uint32_t a, uint32_t b, uint32_t c, uint32_t d, uint32_t subtract) {
    uint8_t a_r = (a >> 24) & 0xFF, a_g = (a >> 16) & 0xFF, a_b = (a >> 8) & 0xFF, a_a = a & 0xFF;
    uint8_t b_r = (b >> 24) & 0xFF, b_g = (b >> 16) & 0xFF, b_b = (b >> 8) & 0xFF, b_a = b & 0xFF;
    uint8_t c_r = (c >> 24) & 0xFF, c_g = (c >> 16) & 0xFF, c_b = (c >> 8) & 0xFF, c_a = c & 0xFF;
    uint8_t d_r = (d >> 24) & 0xFF, d_g = (d >> 16) & 0xFF, d_b = (d >> 8) & 0xFF, d_a = d & 0xFF;
    uint8_t sub_r = (subtract >> 24) & 0xFF, sub_g = (subtract >> 16) & 0xFF, sub_b = (subtract >> 8) & 0xFF, sub_a = subtract & 0xFF;

    uint16_t sum_r = a_r + b_r + c_r + d_r;
    uint16_t sum_g = a_g + b_g + c_g + d_g;
    uint16_t sum_b = a_b + b_b + c_b + d_b;
    uint16_t sum_a = a_a + b_a + c_a + d_a;

    sum_r /= 2;
    sum_g /= 2;
    sum_b /= 2;
    sum_a /= 2;

    uint16_t final_r = sum_r > sub_r ? sum_r - sub_r : 0;
    uint16_t final_g = sum_g > sub_g ? sum_g - sub_g : 0;
    uint16_t final_b = sum_b > sub_b ? sum_b - sub_b : 0;
    uint16_t final_a = sum_a > sub_a ? sum_a - sub_a : 0;

    final_r = final_r > 255 ? 255 : final_r;
    final_g = final_g > 255 ? 255 : final_g;
    final_b = final_b > 255 ? 255 : final_b;
    final_a = final_a > 255 ? 255 : final_a;

    return (final_r << 24) | (final_g << 16) | (final_b << 8) | final_a;
}



void water_ripple_effect(int mouse_pos_x, int mouse_pos_y){
    rafgl_pixel_rgb_t prev_pixel;
    rafgl_pixel_rgb_t curr_pixel;
    rafgl_pixel_rgb_t tmp_pixel;
    int width = upscaled_sea.width;
    int height = upscaled_sea.height;

    /*uint32_t number = 0;

    number |= 000;
    number |= 100 << 8;
    number |= 100 << 16;
    number |= 100 << 24;

    pixel_at_m(upscaled_sea, mouse_pos_x, mouse_pos_y).rgba = number;*/

    // curr je ripple_tmp, prev je ripple_prev, upscaled_sea je tmp (tako mi je bilo najlakse da ga prikazujem posle ove funkcije bez menjanja koda)

    uint16_t damping_mul = 99; // from 100
    uint16_t damping_div = damping_mul; // from 100

    int x,y;
    //printf("TU SAM");

    /*for (x = 1; x < width - 1; x++) {
        for (y = 1; y < height - 1; y++) {
                printf(" upscaled , prev , tmp : %u %u %u \n", pixel_at_m(upscaled_sea, x, y).rgba, pixel_at_m(ripple_prev, x, y).rgba,pixel_at_m(ripple_tmp, x, y).rgba);
        }
    }*/


    for (x = 1; x < width - 1; x++) {
        for (y = 1; y < height - 1; y++) {
            /*if( abs(y * x - mouse_pos_x*mouse_pos_y) < 10){
                printf("Decimal: %u %u %u %u %u\n", pixel_at_m(upscaled_sea, x, y - 1).rgba , pixel_at_m(upscaled_sea, x, y + 1).rgba , pixel_at_m(upscaled_sea, x - 1, y).rgba , pixel_at_m(upscaled_sea, x + 1, y).rgba, pixel_at_m(ripple_prev, x, y).rgba);
            }*/

            uint32_t sum = safe_rgba_operation(pixel_at_m(upscaled_sea, x, y - 1).rgba , pixel_at_m(upscaled_sea, x, y + 1).rgba , pixel_at_m(upscaled_sea, x - 1, y).rgba , pixel_at_m(upscaled_sea, x + 1, y).rgba, pixel_at_m(ripple_prev, x, y).rgba);

            /*if(abs(y * x - mouse_pos_x*mouse_pos_y) < 10){
                printf("Decimal: %u \n", sum);
            }*/
            //uint8_t a_r = (sum >> 24) & 0xFF, a_g = (sum >> 16) & 0xFF, a_b = (sum >> 8) & 0xFF, a_a = sum & 0xFF;
            uint8_t a_r = (sum >> 24) & 0xFF;
            uint8_t a_g = (sum >> 16) & 0xFF;
            uint8_t a_b = (sum >> 8) & 0xFF;
            uint8_t a_a = sum & 0xFF;

            /*if(abs(y * x - mouse_pos_x*mouse_pos_y) < 10){
                printf("CRAZY HAMBURGER:  %u, %u %u %u \n", a_r, a_g, a_b, a_a);
            }*/

            if(abs(y * x - mouse_pos_x*mouse_pos_y) < 10){
                damping_div = 100;
            }else{
                damping_div = damping_mul;
            }
            uint16_t new_r = (a_r * damping_mul) / damping_div;
            uint16_t new_g = (a_g * damping_mul) / damping_div;
            uint16_t new_b = (a_b * damping_mul) / damping_div;
            uint16_t new_a = (a_a * damping_mul) / damping_div;


            /*if(abs(y * x - mouse_pos_x*mouse_pos_y) < 10){
                printf(" NEW:  %u, %u %u %u \n", a_r / damping_mul, a_g * damping_mul , a_b * damping_mul, a_a * damping_mul);
            }*/

            /*if(abs(y * x - mouse_pos_x*mouse_pos_y) < 10){
                printf(" NEW:  %u, %u %u %u \n", new_r, new_g, new_b, new_a);
            }*/
            new_r = new_r > 255 ? 255 : new_r;
            new_g = new_g > 255 ? 255 : new_g;
            new_b = new_b > 255 ? 255 : new_b;
            new_a = new_a > 255 ? 255 : new_a;

            uint32_t final_sum = (new_r << 24) | (new_g << 16) | (new_b << 8) | new_a;

            pixel_at_m(ripple_prev, x , y).rgba = final_sum;
            /*if(abs(y * x - mouse_pos_x*mouse_pos_y) < 10){
                printf("Decimal: %u \n", final_sum);
            }*/

        }
    }


    for (x = 0; x < width; x++) {
        for (y = 0; y < height ; y++) {
            pixel_at_m(ripple_tmp, x, y).rgba = pixel_at_m(upscaled_sea, x, y).rgba;
            pixel_at_m(upscaled_sea, x, y).rgba = pixel_at_m(ripple_prev, x, y).rgba;
            pixel_at_m(ripple_prev, x, y).rgba = pixel_at_m(ripple_tmp, x, y).rgba;
        }
    }
}




typedef struct {
    int x, y;
    int width, height;
} BoundingBox;

int checkCollision(BoundingBox hero, BoundingBox tree) {
    if (hero.x < tree.x + tree.width &&
        hero.x + hero.width > tree.x &&
        hero.y < tree.y + tree.height &&
        hero.y + hero.height > tree.y) {
        return 1;
    }
    return 0;
}

void directional_blur(int hero_pos_x, int hero_pos_y, int direction){
    int x,y;
    float yn,xn;
    int rc,gc,bc;
    int sample_count = 20;
    int vecx = 0, vecy = 0;
    int raster_width = RASTER_WIDTH;
    int raster_height = RASTER_HEIGHT;
    int i,j;
    rafgl_pixel_rgb_t sampled;
    rafgl_pixel_rgb_t result;


    switch(direction) {
        case 0: // Going Down
            vecy = -5;
            vecx = 5;
            break;
        case 1: // Going LEft
            vecx = 20;
            vecy = 5;
            break;
        case 2: // Going Up
            vecy = 20;
            vecx = 5;
            break;
        case 3: // Going Riht
            vecx = -20;
            vecy = 5;
            break;
    }

    for(y = hero_pos_y; y < hero_pos_y + TILE_SIZE; y++)
    {
        yn = 1.0f*y/(hero_pos_y + TILE_SIZE);
        for(x = hero_pos_x; x < hero_pos_x + TILE_SIZE; x++)
        {
            xn = 1.0f*x/(hero_pos_x + TILE_SIZE);
            rc=gc=bc=0;
            for(i = 0; i<sample_count;i++){
                //cist direkcioni blur
                //broj koji polako ide od 0 ka jedinici
                float scalar = 1.0f * i/sample_count;

                //poziciju piksela biramo kao razlicite pozicije na jednom pravcu vektora
                //logika na dalje je ista kao i u prethodnom primeru
                int offsetx = vecx * scalar;
                int offsety = vecy * scalar;

                sampled = pixel_at_m(raster, rafgl_clampi(x+offsetx, 0 , raster_width-1), rafgl_clampi(y+offsety, 0 , raster_height-1));
                rc +=sampled.r;
                gc +=sampled.g;
                bc +=sampled.b;
            }
            rc/=sample_count;
            gc/=sample_count;
            bc/=sample_count;

            result.rgba = rafgl_RGB(rafgl_saturatei(rc),rafgl_saturatei(gc),rafgl_saturatei(bc));
            pixel_at_m(raster,x,y) = result;
        }


    }
}

void flip_hero_and_grayscale(int pos_x_hero, int pos_y_hero, int direction){
    rafgl_pixel_rgb_t sampled;
    rafgl_pixel_rgb_t recolor;
    rafgl_pixel_rgb_t safekeep[TILE_SIZE][TILE_SIZE] = {0};
    int x,y;
        // Create a matrix of pixels
    for (x = pos_x_hero; x < pos_x_hero + TILE_SIZE; x++) {
        for (y = pos_y_hero; y < pos_y_hero + TILE_SIZE; y++) {
            //printf("%d %d ", x,y);
            safekeep[x - pos_x_hero][y - pos_y_hero] = pixel_at_m(raster, x, y);
            //printf("  %d %d \n", x - pos_x_hero, y - pos_y_hero);
        }
        //printf("OUTER LOOP %d \n", y);
    }
    //printf("ANYONE HERE??");

    printf("START  x,y : %d %d  END x,y : %d %d", pos_x_hero, pos_y_hero, pos_x_hero + TILE_SIZE , pos_y_hero +TILE_SIZE);
    for (y = pos_y_hero; y < pos_y_hero + TILE_SIZE; y++) {
        for (x = pos_x_hero; x < pos_x_hero + TILE_SIZE; x++) {

            int sampled_x, sampled_y;

            // Rotate based on direction
            if (direction == 1) {
                // Rotate 90 degrees to the right
                sampled_x = y;
                sampled_y = pos_x_hero + (pos_x_hero + TILE_SIZE - x) ;
                //printf("Hero pos x,y : %d %d  sampled x,y : %d %d , TILE_SIZE, x : %d %d \n", pos_x_hero, pos_y_hero, sampled_x, sampled_y, TILE_SIZE, x);
            } else if (direction == 2) {
                // Rotate 180 degrees
                sampled_x = x;
                sampled_y = pos_y_hero + TILE_SIZE - y + pos_y_hero;
            } else if (direction == 3) {
                // Rotate 90 degrees to the left
                sampled_x = y;
                sampled_y = x;
            } else {
                // No rotation for other directions
                sampled_x = x;
                sampled_y = y;
            }

            if(direction % 2){
                sampled = safekeep[sampled_x - pos_y_hero][sampled_y - pos_x_hero];
            }else{
                sampled = safekeep[sampled_x - pos_x_hero][sampled_y - pos_y_hero];
            }
            //sampled = pixel_at_m(raster,sampled_x,sampled_y);
            int brightness = (sampled.r + sampled.g + sampled.b) / 3;
            //int brightness = (hero.sheet.data.r + hero.sheet.data.g + hero.sheet.data.b) / 3;
            recolor.r = brightness;
            recolor.g = brightness;
            recolor.b = brightness;
            pixel_at_m(raster, x , y) = recolor;
            //printf("%d %d Hero pos x,y : %d %d  sampled x,y : %d %d , camera x,y : %d %d \n",sampled_x - pos_x_hero,sampled_y - pos_y_hero, pos_x_hero, pos_y_hero, sampled_x, sampled_y, camx, camy);
        }
    }
}

void animate_waves(float time){
    int y1, x1, x, y;
    rafgl_pixel_rgb_t sampled;
    int raster_width = upscaled_sea.width;
    int raster_height = upscaled_sea.height;

    //2Pi mi je perioda cosf i sinf sto znaci da time + y treba da ide u rangeu od 2Pi/0.002

    for (y = 0; y < raster_height; y++) {
        for (x = 0; x < raster_width; x++) {
            // Use time directly in the wave motion
            x1 = x + cosf((y + time) * 0.055) * 10; // 0.01 za velike talase (mala frekvencija)
            y1 = y + sinf((x + time) * 0.055) * 10; // 0.01 ...

            sampled = pixel_at_m(upscaled_sea, rafgl_clampi(x1, 0, raster_width - 1), rafgl_clampi(y1, 0, raster_height - 1));
            pixel_at_m(raster, x, y) = sampled;
        }
    }
}

float manual_powf(float base, float exponent) {
    float result = 1.0;

    if (exponent > 0) {
        for (int i = 0; i < exponent; ++i) {
            result *= base;
        }
    } else if (exponent < 0) {
        for (int i = 0; i > exponent; --i) {
            result /= base;
        }
    }

    return result;
}

void fisheye_entrance_effect(int fisheye_cave_x, int fisheye_cave_y){

    int zoom_in_k_frames = number_of_frames;

    float alpha, r;
    float rn;
    //koordinatni pocetak smestamo u centar slike tako sto cemo cx i cy vrednost da oduzimamo i dodajemo na nase koordinate
    int cx=fisheye_cave_x + TILE_SIZE / 2, cy= fisheye_cave_y + TILE_SIZE / 2;
    //poludijagonala(udaljenost nekog coska od centra slike, odnosno koordinatnog pocetka)
    //ona ce nam posluziti za normalizaciju poluprecnika, jer je zapravo i njegova najveca vrednost
    //int diagonal = sqrt((raster_width/2 * raster_width/2) + (raster_height/2 * raster_height/2));

    rafgl_pixel_rgb_t sampled;
    rafgl_pixel_rgb_t result;

    int max_distance_to_corner = 0;

    int top_left = (int)rafgl_distance2D(cx, cy, 0, 0);
    int top_right = (int)rafgl_distance2D(cx, cy, raster_width, 0);
    int bottom_left = (int)rafgl_distance2D(cx, cy, 0, raster_height);
    int bottom_right = (int)rafgl_distance2D(cx, cy, raster_width, raster_height);


    max_distance_to_corner = top_left > top_right ? top_left : top_right;
    max_distance_to_corner = max_distance_to_corner > bottom_left ? max_distance_to_corner : bottom_left;
    max_distance_to_corner = max_distance_to_corner > bottom_right ? max_distance_to_corner : bottom_right;


    int diagonal = max_distance_to_corner;
    //printf("DIJAGONALAAAAAA %d ", diagonal);
    int x,y,y1,x1;

    float inverted_number_of_frames = (float)fisheye_number_of_frames - (float)zoom_in_k_frames;

    float modified_power = 1.0 + 1.0 * inverted_number_of_frames/10;

    //printf(" ZOOM IN K FRAMES, number of frames, inverted_no, modified_pow : %d %d %f %f\n " , zoom_in_k_frames,number_of_frames, inverted_number_of_frames, modified_power);

    printf("%u %u %u", pixel_at_m(raster, 590,416).r,  pixel_at_m(raster, 590,416).g, pixel_at_m(raster, 583,416).b);
    for(y = 0; y < raster_height; y++)
    {
        for(x = 0; x < raster_width; x++)
        {
            ///prebacujemo dekartov sistem koordinata u polarni sistem koordinata
            //funkcija atan2() je zapravo samo funkcija koja radi atan i umesto nas pazi da nemamo deljenje sa 0
            //oduzimamo centar, da bismo vestacki napravili da nam je koordinatni pocetak u centru slike, a ne u gornjem levom uglu
            alpha = atan2(y-cy,x-cx);
            //distanca izmedju centra i nase tacke
            r = rafgl_distance2D(cx,cy, x,y);
            //normalizujemo r
            rn = r/diagonal;

            ///deformacija poluprecnika
            //sa korenovanjem dobijamo udubljeni fisheye
            //rn = sqrtf(rn);
            //sa stepenovanjem dobijamo obrnuti fisheye u odnosu na sqrtf

            //printf("POWF: %f %.3f  %f\n", inverted_number_of_frames/5, powf(rn,1.0 + inverted_number_of_frames/5), rn);

            rn = powf(rn,modified_power); // neki kao fisheye zoom napravljen od broja frejmova i sve vreme se povecava


            ///prebacujemo polarni sistem koordinata u dekartov sistem koordinata
            //preracunavamo r u odnosu na rn koje smo deformisali
            r = rn*diagonal;
            //printf(" cosalfa *r i sinalfa*r : %.3f %3.f\n", cosf(alpha)*r, sin(alpha) * r);
            //nove koordinate x i y, nakon sto smo promenili velicinu r-a
            //potrebno je da ponovno gledamo koordinate u odnosu na pocetak rastera, tj. gornji levi ugao, pa cemo jos i to nadodati
            x1 = cx + cosf(alpha)*r;
            y1 = cy + sinf(alpha)*r;

            if(x==700 && y==500){
                printf("r, rn , diagonal : %f %f %d , cos,sin,alfa: %f %f %f cx,cy %d %d\n",r, rn, diagonal, cos(alpha),sin(alpha), alpha, cx,cy);
                printf("x1,y1 : %d %d\n", x1,y1);
            }



            //semplujemo sa vrednostima novih koordinata clamopovanim u granicama slike
            rafgl_pixel_rgb_t test_pixel = pixel_at_m(raster, x1, y1);

            /*if(x==700 && y==500){
                printf("RGB OF PIXEL %u %u %u\n", test_pixel.r, test_pixel.g, test_pixel.b);
            }*/
            sampled = pixel_at_m(raster, rafgl_clampi(x1, 0, raster_width - 1), rafgl_clampi(y1, 0, raster_height - 1));
            /*if(x >= 700 && x < 710 && y >= 500 && y < 510){
                printf("x1, y1 : %d %d, x,y: %d %d\n", x1,y1, x,y);
                printf("Sampled after rgb: %u %u %u, sampled before r,g,b : %u %u %u\n", sampled.r, sampled.g, sampled.b,test_pixel.r, test_pixel.g, test_pixel.b );
            }*/
            /*if(x==700 && y==500){
                printf("RGB OF PIXEL %u %u %u\n", test_pixel.r, test_pixel.g, test_pixel.b);
                printf("x1,y1 : %d %d RGB OF PIXEL %u %u %u\n", x1,y1 , sampled.r, sampled.g, sampled.b);
            }*/
            pixel_at_m(tmpFisheye, x, y) = sampled;
            /*if(x == 700 && y == 500){
                printf("Number of frames is: %d  %d\n", zoom_in_k_frames, (int)modified_power);
                printf("Which pixel we are looking at  x,y :%d %d x1, y1: %d %d, cam : %d %d\n", x, y, x1, y1, camx, camy);
                printf("cave x, y: %d %d\n ", fisheye_cave_x, fisheye_cave_y);
                printf("cx,cy : %d %d, cos, sin : %.3f %.3f, r, rn , diagonal: %f %f %d\n ", cx, cy, cosf(alpha), sinf(alpha), r, rn , diagonal);
            }*/
        }

    }

     for(y = 0; y < raster_height; y++)
    {
        for(x = 0; x < raster_width; x++)
        {
            pixel_at_m(raster,x,y) = pixel_at_m(tmpFisheye,x,y);
        }
    }
}

void let_there_be_light(int light_source_x[3], int light_source_y[3]){


    int x, y;
    int w = raster_width;
    int h = raster_height;

    rafgl_pixel_rgb_t sampled;
    rafgl_pixel_rgb_t result;
    //poluprecnik kruga
    int r = 350;

    //distanca piksela od centra
    //i koeficijent za jacinu efekta


    for(y = 0; y < h; y++)
    {
        for(x = 0; x < w; x++)
        {
            result = pixel_at_m(raster, x, y);

            int light_idx;

            int cx,cy;
            int collective_sampled_r = 0;
            int collective_sampled_g = 0;
            int collective_sampled_b = 0;
            for (light_idx = 0; light_idx < 3; light_idx++) {
                //printf("%d %d\n", light_source_x[light_idx], light_source_y[light_idx]);
                cx = light_source_x[light_idx];
                cy = light_source_y[light_idx];
                float dist, vignette_factor = 1.5;

                sampled = pixel_at_m(raster, x, y);

                dist = rafgl_distance2D(cx, cy, x, y) / r;
                dist = powf(dist, 2.0f);


                collective_sampled_r += rafgl_clampi(rafgl_saturatei(sampled.r * (1.0f - dist * vignette_factor)), 0 , 255);
                collective_sampled_g += rafgl_clampi(rafgl_saturatei(sampled.g * (1.0f - dist * vignette_factor)), 0 , 255);
                collective_sampled_b += rafgl_clampi(rafgl_saturatei(sampled.b * (1.0f - dist * vignette_factor)), 0 , 255);

            }
            result.r = rafgl_clampi(collective_sampled_r, 0 ,255) / 3;
            result.g = rafgl_clampi(collective_sampled_g, 0 ,255) / 3;
            result.b = rafgl_clampi(collective_sampled_b, 0 ,255) / 3;

            pixel_at_m(raster, x, y) = result;
        }
    }

}


//poziva se samo jednom
void init_tilemap(void)
{
    int x, y;

    //prolazimo kroz svet i randomujemo
    //i imamo 20% ssanse za generisanje drveca i 80%sanse da generisemo travu
    //
    for(y = 0; y < WORLD_SIZE; y++)
    {
        for(x = 0; x < WORLD_SIZE; x++)
        {
            float random_number = 1.0 * rand();
            if(random_number < 0.007 * RAND_MAX)
            {
                tile_world[y][x] = 18;
            }else if(random_number < 0.02 * RAND_MAX)
            {
                tile_world[y][x] = 19;
            }else if (random_number < 0.04 * RAND_MAX){
                tile_world[y][x] = 6;
            }
            else if(random_number < 0.26 * RAND_MAX)
            {
                //generisanje drveca, tajlovi sa indeksima 3,4,5
                tile_world[y][x] = 3 + rand() % 3;
            }
            else
            {
                //generisanje trave, tajlovi sa indeksima 0,1,2
                tile_world[y][x] = rand() % 3;
            }
        }
    }
}

//renderujemo tilemap tako sto ga ucrtavamo u neki raster
//zelimo da iscrtamo tajlove koje samo trenutno mozemo da vidimo a ne sve
//tako da cemo izracunati x0 i y0 kao indeks skroz gornjeg levog tajla
//i x1 i y1 kao indeks skroz donjeg desnog tajla
//zavise nam od pozicije kamere

void render_tilemap(rafgl_raster_t *raster)
{
    int x, y;
    int x0 = camx / TILE_SIZE;
    //koliko nam tajlova staje na ekran tajlovano po x koordinati
    //tj ako nam je ekran sirok 640 piksela, u njega staje 10 tajlova, a ovo +1 dodajemo
    //u slucaju kada se nalazimo na nekoj ivici, pa nam se prikazuje po pola tajla
    //i za y imamo +2, jer ima onih drveca koji imaju viska gore
    //pa da ne bismo to odsekli
    int x1 = x0 + (raster_width / TILE_SIZE) + 1;
    int y0 = camy / TILE_SIZE;
    int y1 = y0 + (raster_height / TILE_SIZE) + 2;

    //klempujemo
    if(x0 < 0) x0 = 0;
    if(y0 < 0) y0 = 0;
    if(x1 < 0) x1 = 0;
    if(y1 < 0) y1 = 0;

    //klempujemo
    if(x0 >= WORLD_SIZE) x0 = WORLD_SIZE - 1;
    if(y0 >= WORLD_SIZE) y0 = WORLD_SIZE - 1;
    if(x1 >= WORLD_SIZE) x1 = WORLD_SIZE - 1;
    if(y1 >= WORLD_SIZE) y1 = WORLD_SIZE - 1;

    //pointer ka tipu tajla kojeg trenutno iscrtavamo
    //inace tajlovi su nam ucitani dole u for petlji u inicijalizaciji
    rafgl_raster_t *draw_tile;

    //prolazimo sa y i x od y0 do y1 i od x0 do x1, to su svi tajlovi koji se vide
    //tj svi tajlovi koje treba da iscrtamo
    //kako izvlacimo koji je to tip tajla, tako sto uzimamo tile_world sa tom pozicijom i to modujemo sa brojem tajlova
    //i zatim to dodoao na niz tiles, tj pokazivac na pocetak niza
    // ovo je kao da smo uradili &tiles[tile_world[y][x] % NUMBER_OF_TILES)]
    //zatim pozivamo funkciju za iscrtavanje rastera i njoj prosledjujemo na koji raster, kog tipa je tajl i
    //na kojoj poziciji
    //x*tile_size(tj to je pozicija tog tajla u svetu), ali mi ga zelimo u odnosu na pocetak ekrana, tako da zato oduzimamo i camx
    //
    //u vecini slucajeva ce nam ovo (- draw_tile->height + TILE_SIZE) biti 0, tj to nece raditi nista
    // ali u slucaju da je taj tajl malo veci po y, njega cemo samo pomeriti malo gore
    //
    //nama je bitan redosled iscrtavanja tajlova, tako da moramo prvo da iscrtamo sve gornje, pa tek onda sve donje tajlove
    //jer donji moze da sadrzi ono sto se crta preko onog iznad
    //da radimo suprotno, ne bismo videli vrhove jelki, jer bi one bile ispod tajlova koji su se preslikali preko njih
    //
    //a ukoliko kod y ne bismo stavili ovaj offset, tajl bi krenuo da se iscrtava odakle bi trebao inace, tako da bismo imali vrhove i providan deo
    //dok bi donji deo jelkica bio presecen, jer se nista ne bi iscrtavalo nakon sto se iscrta prvih 64piksela po y
    //
    for(y = y0; y <= y1; y++)
    {
        for(x = x0; x <= x1; x++)
        {
            draw_tile = tiles + (tile_world[y][x] % NUMBER_OF_TILES);
            rafgl_raster_draw_raster(raster, draw_tile, x * TILE_SIZE - camx, y * TILE_SIZE - camy - draw_tile->height + TILE_SIZE);
        }
    }
    if(fisheye_entrance && number_of_frames > 0){
        fisheye_entrance_effect(fisheye_cave_x - camx,fisheye_cave_y - camy);
        number_of_frames--;

        //sleep(10);//printf(" U IFUUUUU SMOOOO : entrace , no frames, fisheye_cave_x - camx, fisheye_cave_y - camy %d %d %d %d\n" , fisheye_entrance, number_of_frames, fisheye_cave_x - camx ,fisheye_cave_y - camy);
    }

    rafgl_raster_draw_rectangle(raster, selected_x * TILE_SIZE - camx, selected_y * TILE_SIZE - camy, TILE_SIZE, TILE_SIZE, rafgl_RGB(255, 255, 0));

}


/*int detect_collision(int x, int y){
    //printf("x, y : %d, %d", x,y);
    int first_coord = rafgl_clampi((x)/ TILE_SIZE, 0, WORLD_SIZE - 1);
    int second_coord = rafgl_clampi((y) / TILE_SIZE, 0, WORLD_SIZE - 1);

    //printf("First, Second: %d, %d\n", first_coord,second_coord);
    if(tile_world[first_coord][second_coord] < 3){
        //printf("NIJE POPUNJENA %d %d : %d \n", first_coord,second_coord,tile_world[first_coord][second_coord] );
        return 0;
    }
    printf("JESTE POPUNJENA SA DRVETOM %d, %d  : %d \n", x, y, tile_world[first_coord][second_coord] );
    return 1;
}*/




void main_state_init(GLFWwindow *window, void *args, int width, int height)
{

    //printf("IZ MAINA x,y  : %d  %d \n", hero_pos_x, hero_pos_y );

    //int j = detect_collision(hero_pos_x,hero_pos_y);
    //while(j){
       // hero_pos_x += TILE_SIZE;
       // camx += TILE_SIZE;
        //j = detect_collision(hero_pos_x,hero_pos_y);
    //}
    //printf("IZ MAINA x,y  : %d  %d \n", hero_pos_x, hero_pos_y );
    /* inicijalizacija */
    /* raster init nam nije potreban ako radimo load from image */

    time = 1;
    ripple_initialized = 0;
    frames = 0;
    is_cave = 0;
    fisheye_entrance = 0;
    turned_on_lights = 0;
    hero_on_ground = 1;
    last_clicked = 0;
    hero_vertical_velocity = 0;
    ground_level = hero_pos_cave_y;


    //kamera x i kamera y je postavljeno na 600,600 na pocetku
    int camx = 0, camy = 0;

    int hero_is_dead = 0;
    rafgl_raster_load_from_image(&sea, "res/images/talasi1.png");
    rafgl_raster_load_from_image(&checker, "res/images/checker32.png");

    rafgl_raster_load_from_image(&cave_instance1, "res/images/inside_cave1.png");
    rafgl_raster_load_from_image(&cave_instance2, "res/images/inside_cave2.png");

    raster_width = width;
    raster_height = height;

    rafgl_raster_init(&upscaled_sea, raster_width, raster_height);
    rafgl_raster_bilinear_upsample(&upscaled_sea, &sea);

    rafgl_raster_init(&ripple_prev, raster_width, raster_height);
    rafgl_raster_bilinear_upsample(&ripple_prev, &sea);

    rafgl_raster_init(&ripple_tmp, raster_width, raster_height);
    rafgl_raster_bilinear_upsample(&ripple_tmp, &sea);


    rafgl_raster_init(&cave_instance1_scaled, raster_width, raster_height);
    rafgl_raster_bilinear_upsample(&cave_instance1_scaled, &cave_instance1);

    rafgl_raster_init(&cave_instance2_scaled, raster_width, raster_height);
    rafgl_raster_bilinear_upsample(&cave_instance2_scaled, &cave_instance2);

    rafgl_raster_init(&tmpFisheye, raster_width, raster_height);
    rafgl_raster_bilinear_upsample(&tmpFisheye, &sea);


    rafgl_raster_init(&raster, raster_width, raster_height);
    rafgl_raster_init(&raster2, raster_width, raster_height);


    int i;


    //naziv putanje do tajla kog zelimo da ucitamo
    char tile_path[256];

    //prolazimo kroz sve tajlove i ucitavamo slike u rastere za odgovarajuce pozicije u nizu
    for(i = 0; i < NUMBER_OF_TILES; i++)
    {
        //ovde ce see samo zameniti ovo %d sa trenutnim indeksom u foru
        //tako da cemo u tiles niz imati ucitane sve slike redosledom od 0 do16
        sprintf(tile_path, "res/tiles/svgset%d.png", i);
        rafgl_raster_load_from_image(&tiles[i], tile_path);
    }

    //spritesheet je inicijalizovan tako sto prosledimo
    //koju sliku ucitavamo
    //i koliko ima slicica po x, a koliko po y
    //i to smo samo zapamtili u sprite sheet hero
    rafgl_spritesheet_init(&hero, "res/images/character.png", 10, 4);

    rafgl_spritesheet_init(&hero1, "res/images/character1.png", 8, 4);

    //i sada kada imamo ucitane tajloce mozemo i da inicijalizujemo tilemapu
    init_tilemap();
    /*for (int i = 0; i < WORLD_SIZE; i++) {
        for (int j = 0; j < WORLD_SIZE; j++) {
            printf("%d", tile_world[i][j]);
        }
        printf("\n");
    }*/
    fisheye_cave_x = RASTER_WIDTH/2;
    fisheye_cave_y = RASTER_HEIGHT/2;
    number_of_frames = 30;
    fisheye_number_of_frames = number_of_frames;

    rafgl_texture_init(&texture);
    int x,y;
    for(y = 0 ; y < WORLD_SIZE; y++){
        for(x = 0 ; x  < WORLD_SIZE; x++){
            printf("%d ", tile_world[y][x]);
        }
        printf("\n");
    }
}


int pressed;
float location = 0;
float selector = 0;

//flag da li se krecemo ili stojimo
int animation_running = 0;
//koji je indeks slike po x osi iz spritesheeta, tj koju kolonu gledamo
int animation_frame = 0;
//koji je indeks slike po y osi iz spritesheeta, tj koji red gledamo
int direction = 0;


//brzina kretanja karaktera
int hero_speed = 150;

int flipped = 0;

//da se ne menja u svakom frameu, da ne bismo prebrzo isli kroz animaciju
//vec kazemo koliko frameova da traje neka animacija, tj jedna slicica
int hover_frames = 0;

void main_state_update(GLFWwindow *window, float delta_time, rafgl_game_data_t *game_data, void *args)
{
    /* hendluj input */
    if(game_data->is_lmb_down && game_data->is_rmb_down)
    {
        pressed = 1;
        location = rafgl_clampf(game_data->mouse_pos_y, 0, raster_height - 1);
        selector = 1.0f * location / raster_height;
    }
    else
    {
        pressed = 0;
    }

    if(game_data->keys_down[RAFGL_KEY_LEFT_SHIFT]){
        hero_speed = 250;
    }else{
        hero_speed = 150;
    }


    //racunanje u kom polju mape se nalazi mis
    //racunanje je objasnjeno na pocetku koda
    //klempovano izmedju 0 i world_sizea
    //kada ne bismo imali ovo klempovanje, samo bismo mogli da selektujemo i van worlda
    selected_x = rafgl_clampi((game_data->mouse_pos_x + camx) / TILE_SIZE, 0, WORLD_SIZE - 1);
    selected_y = rafgl_clampi((game_data->mouse_pos_y + camy) / TILE_SIZE, 0, WORLD_SIZE - 1);

    /* izmeni raster */

       //pretpostavimo da se krecemo
    animation_running = 1;

    int tmp_hero_pos_x = hero_pos_x;
    int tmp_hero_pos_y = hero_pos_y;

    int delta_pos_x = 0;
    int delta_pos_y = 0;

    //printf("ENTRANCE AND NUMBER_OF_FRAMES : %d %d" , fisheye_entrance, number_of_frames);
    /*if(fisheye_entrance && number_of_frames > 0){
        fisheye_entrance_effect(fisheye_cave_x - camx,fisheye_cave_y - camy);
        number_of_frames--;
        //printf(" U IFUUUUU SMOOOO : entrace , no frames, fisheye_cave_x - camx, fisheye_cave_y - camy %d %d %d %d\n" , fisheye_entrance, number_of_frames, fisheye_cave_x - camx ,fisheye_cave_y - camy);
    }*/


    if(hero_is_dead){
        if(!flipped){
            //printf("%d %d , Hero pos x,y " , hero_pos_x, hero_pos_y);
            flip_hero_and_grayscale(hero_pos_x - camx ,hero_pos_y - camy,direction);
            flipped = 1;
        }
    }else if(!is_cave || number_of_frames > 0){


        if(fisheye_entrance || !is_cave){
            // smanjujem ga ovde da bih tacno usao 30 puta u elif umesto 29 (u slucaju da sam u gornjem ifu stavio) - ipak kvari logiku xd

            //direkcija nam predstavlja koji red iz spritesheeta treba da uzmemo, tj to nam je kao y iz spitesheeta
            //broj koji dodeljujemo directionu, zavisi od organizacije pravaca, tj redova u okviru spritesheeta

            //ako je duplo w pritisnuto, onda menjamo poziciju heroa tako sto ga pomeramo ka gore i stavimo da je njegova direkcija 2
            if(game_data->keys_down[RAFGL_KEY_W])
            {
                delta_pos_y = -1.0 * (hero_speed * delta_time);
                hero_pos_y = hero_pos_y + delta_pos_y;
                direction = 2;
            }
             //ako je s pritisnuto, onda menjamo poziciju heroa tako sto ga pomeramo ka dole i stavimo da je njegova direkcija 0
            else if(game_data->keys_down[RAFGL_KEY_S])
            {
                delta_pos_y = hero_speed * delta_time;
                hero_pos_y = hero_pos_y + delta_pos_y;
                direction = 0;
            }
            //ako je s pritisnuto, onda menjamo poziciju heroa tako sto ga pomeramo ka levo i stavimo da je njegova direkcija 1
            else if(game_data->keys_down[RAFGL_KEY_A])
            {
                delta_pos_x = -1.0 * (hero_speed * delta_time);
                hero_pos_x = hero_pos_x + delta_pos_x;
                direction = 1;
            }
            //ako je s pritisnuto, onda menjamo poziciju heroa tako sto ga pomeramo ka desno i stavimo da je njegova direkcija 3
            else if(game_data->keys_down[RAFGL_KEY_D])
            {
                delta_pos_x = hero_speed * delta_time;
                hero_pos_x = hero_pos_x + delta_pos_x;
                direction = 3;
            }
            else
            {
                //ipak se ne krecemo
                animation_running = 0;
            }

            hero_pos_x = rafgl_clampi(hero_pos_x, 0, (WORLD_SIZE - 1) * TILE_SIZE);
            hero_pos_y = rafgl_clampi(hero_pos_y, 0, (WORLD_SIZE - 1) * TILE_SIZE);

            int does_it_collide = 0;
            int x,y;
            BoundingBox heroz = {hero_pos_x,hero_pos_y,TILE_SIZE - 30,TILE_SIZE - 30}; // smanjujem malo heroja da moze da prodje (da se ne klipuje)
            //krecem od gornjeg desnog ugla i velicina mu je 64-20, 64-20, tako i racunam collisione
            //printf("%d , %d : Hero x,y " , hero_pos_x, hero_pos_y);

            //testiranje za svih 8 susednih piksela (u realnosti treba 4 ali nema veze, nisam primetio na pocektu da nema dijagonalnog kretanja)
            for ( x = - 1; x <= 1; x++ ) {
                for (y = - 1; y <= 1; y++) {
                    int first_coord = rafgl_clampi(hero_pos_y / TILE_SIZE + x, 0, WORLD_SIZE - 1);
                    int second_coord = rafgl_clampi(hero_pos_x / TILE_SIZE + y, 0, WORLD_SIZE - 1);
                    int curr_x = rafgl_clampi(hero_pos_x / TILE_SIZE, 0, WORLD_SIZE - 1);
                    int curr_y = rafgl_clampi(hero_pos_y / TILE_SIZE, 0, WORLD_SIZE - 1);

                    //printf("curr : %d , %d , FIRST second  coord : %d, %d  ----  tile_World val: %d \n",  curr_x, curr_y, second_coord, first_coord, tile_world[first_coord][second_coord]);


                    //printf("CHECK0 CHECK0\n");

                    /*for (int i = 0; i < WORLD_SIZE; i++) {
                        for (int j = 0; j < WORLD_SIZE; j++) {
                        printf("%d", tile_world[i][j]);
                        }
                        printf("\n");
                    }*/


                    //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                    if (tile_world[first_coord][second_coord] >= 3 && tile_world[first_coord][second_coord] <= 6 || tile_world[first_coord][second_coord] == 18 || tile_world[first_coord][second_coord] == 19) {
                        BoundingBox tree = {second_coord * TILE_SIZE, first_coord * TILE_SIZE, TILE_SIZE, TILE_SIZE};
                        //printf("STAAAAAA %d %d\n", second_coord * TILE_SIZE, first_coord * TILE_SIZE);
                        //printf("CHECK1 CHECK1\n");
                        if (checkCollision(heroz, tree)) {
                            //printf("CHECK2 CHECK2\n");
                            does_it_collide = 1;
                            if( tile_world[first_coord][second_coord] == 6){
                                hero_is_dead = 1;
                            }else if( tile_world[first_coord][second_coord] == 18 || tile_world[first_coord][second_coord] == 19){
                                fisheye_entrance = 1;
                                is_cave = 1;
                                fisheye_cave_x = second_coord * TILE_SIZE;
                                fisheye_cave_y = first_coord * TILE_SIZE;

                                //printf("MY LALALAAJF %d %d\n", fisheye_cave_x, fisheye_cave_y);
                                //fisheye_entrance_effect(fisheye_cave_x - camx , fisheye_cave_y - camy);
                                //printf("TUUUaaa SMO");
                            }
                        }
                    }
                }
            }

            if(!does_it_collide){
                camx += delta_pos_x;
                camy += delta_pos_y;
            }else{
                hero_pos_x = tmp_hero_pos_x;
                hero_pos_y = tmp_hero_pos_y;
            }

            //printf("tmp x, y: %d, %d \n", tmp_hero_pos_x, tmp_hero_pos_y);
            //printf("Selected x, y : %d, %d  ----  cam x, y: %d, %d \n", hero_pos_x, hero_pos_y, camx, camy);


            //ako zelimo da animiramo, tj kliknuto je nesto od wasd, dakle krecemo se, pa zato i hocemo da promenimo polozaj
            if(animation_running)
            {
                //hover frames je koliko frameova zelimo da stojimo na jednom spriteu, tj svaki sprite ce se zadrzati po 5 frameova
                //i tek kad hover padne na 0 mi cemo preci na sledecu animaciju
                if(hover_frames == 0)
                {
                    //i samo kazemo ovo, tj kazemo da se kroz ove animacije krecemo na desno
                    //i ako smo dosli do kraja, da se vratimo na pocetku tog reda
                    //tj animation frame nam je x pozicija u jednom redu dokle smo stigli, kao x iz spritesheeta
                    animation_frame = (animation_frame + 1) % 10;
                    //kada bismo ovde stavili 1, menjali bismo animaciju kretanja u svakom frameu, pa bi nam bilo prebrzo
                    hover_frames = 5;
                }
                else
                {
                    //smanjujemo hover dok ne dodjemo do 0
                    hover_frames--;
                }
            }

            if(game_data ->is_lmb_down){
                frames = 5;
                uint32_t number = 0;

                number |= 255;
                number |= 255 << 8;
                number |= 255 << 16;
                number |= 255 << 24;

                x_pos_mouse = game_data->mouse_pos_x;
                y_pos_mouse = game_data->mouse_pos_y;


               pixel_at_m(upscaled_sea, x_pos_mouse, y_pos_mouse).rgba = number;

               pixel_at_m(upscaled_sea, rafgl_clampi(x_pos_mouse - 5, 0, raster_width),rafgl_clampi(y_pos_mouse - 5, 0, raster_height)).rgba = number;

               pixel_at_m(upscaled_sea, rafgl_clampi(x_pos_mouse + 5, 0, raster_width),rafgl_clampi(y_pos_mouse + 5, 0, raster_height)).rgba = number;
            }


            if(frames > 0){
                water_ripple_effect(game_data->mouse_pos_x, game_data->mouse_pos_y);
                frames--;
            }
        }
        float xn, yn;
        int x,y;

        rafgl_pixel_rgb_t sampled, sampled2, resulting, resulting2;

        for(y = 0; y < raster_height; y++)
        {
            yn = 1.0f * y / raster_height;
            for(x = 0; x < raster_width; x++)
            {
                xn = 1.0f * x / raster_width;

                sampled = pixel_at_m(upscaled_sea, x, y);
                sampled2 = rafgl_point_sample(&sea, xn, yn);

                resulting = sampled;
                resulting2 = sampled2;

                resulting.rgba = rafgl_RGB(0, 0, 0);
                resulting = sampled;

                pixel_at_m(raster, x, y) = resulting;
                pixel_at_m(raster2, x, y) = resulting2;


                if(pressed && rafgl_distance1D(location, y) < 3 && x > raster_width - 15)
                {
                    pixel_at_m(raster, x, y).rgba = rafgl_RGB(255, 0, 0);
                }

            }
        }

        time += 1;
        // objasnjeno zasto gore u funkciji animate waves
        if(time > 3141){
            time = 1;
        }

        animate_waves(time);

        //kada se pritisne plus da se samo menja tajl koji je selektovan
        //menjamo brojeve u matrici tile_world od kojih nam zavisi koji tajl cemo da setujemo na tu poziciju
        if(game_data->keys_pressed[RAFGL_KEY_KP_ADD])
        {
            tile_world[selected_y][selected_x]++;
            tile_world[selected_y][selected_x] %= NUMBER_OF_TILES;
        }


        //na kraju samo pozovemo i renderovanje tog rastera
        render_tilemap(&raster);


        //pomocna funkcija koja kaze na koji raster da ucrta spritesheet
        //a spritesheet ce u initu biti veoma jednostavno definisan
        //kada zelimo da iscrtamo sprite hero, samo cemo staviti u koji raster iscrtavamo
        //koji spritesheet, koji x frame, koji y frame
        //i na kojoj poziciji unutar tog rastera zelimo da iscrtamo
        rafgl_raster_draw_spritesheet(&raster, &hero, animation_frame, direction, hero_pos_x - camx, hero_pos_y - camy);

        if(hero_speed == 250){
            directional_blur(hero_pos_x - camx , hero_pos_y - camy, direction);
        }

        /* shift + s snima raster */
        //za skrinsotovanje ekrana
        if(game_data->keys_pressed[RAFGL_KEY_S] && game_data->keys_down[RAFGL_KEY_LEFT_SHIFT])
        {
            sprintf(save_file, "save%d.png", save_file_no++);
            rafgl_raster_save_to_png(&raster, save_file);
        }
    }else{

        number_of_frames = 0; // ne sme odmah da se updatuje number of frames jer je povezan sa entrancom (uci ce u drugi if)
        // tek kad izlazimo updatujemo ovu promenljivu iz cave-a
        fisheye_entrance = 0; //resetujem ovo da ne bih ulazio u if iznad svega opet

        //printf(" TU SMO ");
        //ako je s pritisnuto, onda menjamo poziciju heroa tako sto ga pomeramo ka levo i stavimo da je njegova direkcija 1
        if(game_data->keys_down[RAFGL_KEY_A])
        {
            hero_pos_cave_x = hero_pos_cave_x - hero_speed * delta_time;
            direction = 2;
        }
        //ako je s pritisnuto, onda menjamo poziciju heroa tako sto ga pomeramo ka desno i stavimo da je njegova direkcija 3
        else if(game_data->keys_down[RAFGL_KEY_D])
        {
            hero_pos_cave_x = hero_pos_cave_x + hero_speed * delta_time;
            direction = 3;
        }else if(game_data->keys_down[RAFGL_KEY_SPACE]){
            if(last_clicked == RAFGL_KEY_SPACE){
                hero_vertical_velocity += JUMP_VELOCITY/8;
                if(hero_pos_cave_x > 0 && hero_pos_cave_x < raster_width/3){
                    turned_on_lights = 1;
                }
            }
            if (hero_on_ground) {
                hero_vertical_velocity = JUMP_VELOCITY;
                hero_on_ground = 0;
                last_clicked = RAFGL_KEY_SPACE;
            }
            direction = 2;
        }
        else
        {
            //ipak se ne krecemo
            animation_running = 0;
        }

        //ako zelimo da animiramo, tj kliknuto je nesto od wasd, dakle krecemo se, pa zato i hocemo da promenimo polozaj
        if(animation_running)
        {
            //hover frames je koliko frameova zelimo da stojimo na jednom spriteu, tj svaki sprite ce se zadrzati po 5 frameova
            //i tek kad hover padne na 0 mi cemo preci na sledecu animaciju
            if(hover_frames == 0)
            {

                //i samo kazemo ovo, tj kazemo da se kroz ove animacije krecemo na desno
                //i ako smo dosli do kraja, da se vratimo na pocetku tog reda
                //tj animation frame nam je x pozicija u jednom redu dokle smo stigli, kao x iz spritesheeta
                animation_frame = (animation_frame + 1) % 5;
                //kada bismo ovde stavili 1, menjali bismo animaciju kretanja u svakom frameu, pa bi nam bilo prebrzo
                /*if(hero_vertical_velocity > 0){
                    anima
                }*/
                hover_frames = 5;
            }
            else
            {
                //smanjujemo hover dok ne dodjemo do 0
                hover_frames--;
            }

        }

        // Apply gravity and update the vertical position
        hero_vertical_velocity -= GRAVITY * delta_time;
        hero_pos_cave_y -= hero_vertical_velocity * delta_time;

        // Check if the hero is on the ground (you may need to adjust the value for "ground_level")
        if (hero_pos_cave_y >= ground_level) {
            hero_pos_cave_y = ground_level; // Snap the hero to the ground
            hero_on_ground = 1; // Set the flag to indicate that the hero is on the ground
            hero_vertical_velocity = 0.0f; // Reset the vertical velocity
        }

        if(hero_pos_cave_y <= 0){
            hero_pos_cave_y = 3;
            hero_vertical_velocity = 0;
        }

        // izmeni raster

        int x, y;

        float xn, yn;

        rafgl_pixel_rgb_t sampled, sampled2, resulting, resulting2;


        for(y = 0; y < raster_height; y++)
        {
            yn = 1.0f * y / raster_height;
            for(x = 0; x < raster_width; x++)
            {
                xn = 1.0f * x / raster_width;

                sampled = pixel_at_m(cave_instance1, x, y);

                resulting = sampled;

                resulting.rgba = rafgl_RGB(0, 0, 0);

                resulting = sampled;

                pixel_at_m(cave_instance1, x, y) = resulting;

            }
        }

        //(raster, draw_tile, x * TILE_SIZE - camx, y * TILE_SIZE - camy - draw_tile->height + TILE_SIZE);
        rafgl_raster_draw_raster(&raster, &cave_instance1_scaled, 0, 0);
        rafgl_raster_draw_spritesheet(&raster, &hero1, animation_frame, direction, hero_pos_cave_x, hero_pos_cave_y);
        if(!turned_on_lights){
            int light_source_x[3]; // nije radilo da u argumentu dam MAX_LIGHT_SOURCES, nice :)
            int light_source_y[3];
            int i = 0;
            for(i = 0; i < 3; i++){
                light_source_x[i] = (int)(raster_width/6) +  i*((int)(5 * raster_width / 12));
                light_source_y[i] = (int)raster_height/6;
                //printf("X Y i LS: %d %d %d %d\n", light_source_x[i], light_source_y[i], i, i*((int)(5 /12 * raster_width)));
            }
            let_there_be_light(light_source_x, light_source_y);
        }
        //pomocna funkcija koja kaze na koji raster da ucrta spritesheet
        //a spritesheet ce u initu biti veoma jednostavno definisan
        //kada zelimo da iscrtamo sprite hero, samo cemo staviti u koji raster iscrtavamo
        //koji spritesheet, koji x frame, koji y frame
        //i na kojoj poziciji unutar tog rastera zelimo da iscrtamo

    }
}


void main_state_render(GLFWwindow *window, void *args)
{
    /* prikazi teksturu */
    rafgl_texture_load_from_raster(&texture, &raster);
    rafgl_texture_show(&texture, 0);
}


void main_state_cleanup(GLFWwindow *window, void *args)
{
    rafgl_raster_cleanup(&raster);
    rafgl_raster_cleanup(&raster2);
    rafgl_raster_cleanup(&upscaled_sea);
    rafgl_raster_cleanup(&ripple_prev);
    rafgl_raster_cleanup(&ripple_tmp);
    rafgl_texture_cleanup(&texture);

}
