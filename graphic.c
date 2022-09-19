#include "dhex_window_app.h"
#include <stdio.h>
#include <stdlib.h>

#define SCREEN_W 200
#define SCREEN_H 200

#define MICROS_PER_FRAME (16667 * 2)

//resolución de los heat_plate(s) en su formato gráfico
#define GRAPHIC_N 200

static inline int iclamp(int v, int vmin, int vmax) {
	if (v < vmin) return vmin;
	if (v > vmax) return vmax;
	return v;
}

typedef struct {
	// access with [x][y]
	uint8_t cell[GRAPHIC_N][GRAPHIC_N];
} graphic_heat_plate;

static void copy_graphic_heat_plate_into_dhex_sprite(graphic_heat_plate* ghp, dhex_sprite* spr) {
	if (spr->w != GRAPHIC_N || spr->h != GRAPHIC_N) return;
	for (int py = 0; py < GRAPHIC_N; py++)
	for (int px = 0; px < GRAPHIC_N; px++) {
		spr->pixeldata[px + GRAPHIC_N * py] = ghp->cell[px][py];
	}
}

//parameter: filename for visualization
int main(int argc, char** argv) {
	if (argc < 2) {
		printf("Usage: graphic.exe <data filename>\n");
		return 0;
	}
	
	FILE* f = fopen(argv[1], "rb");
	if (!f) {
		printf("Couldn't load data file\n");
		return 0;
	}
	
	int16_t nframes;
	fread(&nframes, sizeof(int16_t), 1, f);
	if (nframes <= 0) {
		printf("Bad data file\n");
		return 0;
	}
	
	graphic_heat_plate** graphic_history = malloc(sizeof(graphic_heat_plate*) * nframes );
	
	for (int k = 0; k < nframes; k++) {
		graphic_history[k] = malloc(sizeof(graphic_heat_plate));
		for (int py = 0; py < GRAPHIC_N; py++)
		for (int px = 0; px < GRAPHIC_N; px++) {
			fread(&(graphic_history[k]->cell[px][py]), sizeof(uint8_t), 1, f);
		}
	}
	
	dhex_app* app = dhex_create_app("Heated plate visualizer", SCREEN_W, SCREEN_H, 4.0/3.0, 1, 0);
	dhex_set_window_size(app, SCREEN_W * 2, SCREEN_H * 2);
	
	int current_frame = 0;
	int delta_frame = 0;
	
	dhex_palette pal;
	pal.ncolors = 255;
	for (int k = 0; k <= 255; k++) {
		pal.colors[k] = (dhex_rgba){(uint8_t)k, 0, 0, 255};
	}
	
	dhex_sprite frame;
	frame.w = GRAPHIC_N;
	frame.h = GRAPHIC_N;
	frame.pixeldata = malloc(sizeof(uint8_t) * GRAPHIC_N * GRAPHIC_N);
	
	while (app->quit == 0) {
		dhex_set_timer();
		dhex_refresh_app(app);
		
		if (app->input.mouse_pressed) {
			current_frame = 0;
			delta_frame = 1;
		}
		
		dhex_clear_canvas(app->main_canvas);
		copy_graphic_heat_plate_into_dhex_sprite(graphic_history[current_frame], &frame);
		dhex_draw_sprite(&frame, &pal, app->main_canvas, 0, 0, 0, 0);
		dhex_render_canvas(app->main_canvas); 
		
		current_frame += delta_frame;
		current_frame = iclamp(current_frame, 0, nframes - 1);
		
		dhex_refresh_window(app); 
		dhex_wait_timer(MICROS_PER_FRAME);
	}
	
	dhex_free_app(app);
	return 0;

}