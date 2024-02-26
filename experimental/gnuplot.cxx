// Gnuplot.cxx
// V 1.06
//
// 04/10/22: addition of "\n\n" at "set term jpeg" (reason unknown)

// --------------- Header ----------------------------------------------

#include <string>

typedef std::string string;

void gnuplot_one_function (string title, string style, string label_x, string label_y, double *x, double *y, int n);
void gnuplot_one_function_jpg (string title, string style, string label_x, string label_y, double *x, double *y, int n, string filename);
void gnuplot_one_function_square_jpg (string title, string style, string label_x, string label_y, double *x, double *y, int n, string filename);
void gnuplot_one_function_general (string title, string style, string label_x, string label_y, double *x, double *y, int n, bool jpg_file, string filename, bool square, bool ranges, double xmin, double xmax, double ymin, double ymax);

void gnuplot_one_function_3d (string title, string style, string label_x, string label_y, string label_z, double *x, double *y, void *z, int nx, int ny, int connected, double view1, double view2);
void gnuplot_one_function_3d_jpg (string title, string style, string label_x, string label_y, string label_z, double *x, double *y, void *z, int nx, int ny, int connected, double view1, double view2, string filename);

void gnuplot_two_functions (string title, string style, string label_x, string label_y, double *x, double *y, int n, string dataname, double *x2, double *y2, int n2, string dataname2);
void gnuplot_two_functions_jpg (string title, string style, string label_x, string label_y, double *x, double *y, int n, string dataname, double *x2, double *y2, int n2, string dataname2, string filename);
void gnuplot_two_functions_square_jpg (string title, string style, string label_x, string label_y, double *x, double *y, int n, string dataname, double *x2, double *y2, int n2, string dataname2, string filename);
void gnuplot_two_functions_general (string title, string label_x, string label_y, double *x, double *y, int n, string dataname, string style, double *x2, double *y2, int n2, string dataname2, string style2, bool jpeg, string filename, bool square, bool ranges, double xmin, double xmax, double ymin, double ymax);

void gnuplot_two_arrays_in_datafile (double *x, double *y, int n, string filename);


// ----------------- MAIN ----------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#ifdef WIN32
#include <io.h>
#endif // #ifdef WIN32


FILE *gnuplot_tmp_file;
FILE *gnuplot_command;



string gnuplot_tmp_file_name = "gnuplot_tmp_file.tmp";


void gnuplot_init (void)
{
	gnuplot_command = popen ("gnuplot", "w");
	if (gnuplot_command == NULL)
	{
		fprintf (stderr, "Error starting gnuplot (maybe gnuplot path is not set)!\n");
		exit (EXIT_FAILURE);
	}
}


void gnuplot_open_tmp_file ()
{
	gnuplot_tmp_file = fopen (gnuplot_tmp_file_name.c_str (), "w");
	if (gnuplot_tmp_file == NULL)
	{
		fprintf (stderr, "Error opening tmp gnuplot file!\n");
		exit (EXIT_FAILURE);
	}
}


void gnuplot_close_tmp_file ()
{
	fclose (gnuplot_tmp_file);
}


void gnuplot_remove_tmp_file ()
{
	remove (gnuplot_tmp_file_name.c_str ());
}


void gnuplot_close ()
{
	if (pclose (gnuplot_command) == -1)
	{
		fprintf (stderr, "Closing communication to gnuplot not possible!\n");
		exit (EXIT_FAILURE);
	}
}


void gnuplot_set_title (string label)
{
	fprintf (gnuplot_command, "set title \"%s\"\n", label.c_str ());
}


void gnuplot_set_xlabel (string label)
{
	fprintf (gnuplot_command, "set xlabel \"%s\"\n", label.c_str ());
}


void gnuplot_set_ylabel (string label)
{
	fprintf (gnuplot_command, "set ylabel \"%s\"\n", label.c_str ());
}


void gnuplot_set_zlabel (string label)
{
	fprintf (gnuplot_command, "set zlabel \"%s\"\n", label.c_str ());
}


void gnuplot_set_xrange (double min, double max)
{
	fprintf (gnuplot_command, "set xrange [%.2f:%.2f]\n", min, max);
}


void gnuplot_set_yrange (double min, double max)
{
	fprintf (gnuplot_command, "set yrange [%.2f:%.2f]\n", min, max);
}


void gnuplot_set_zrange (double min, double max)
{
	fprintf (gnuplot_command, "set zrange [%.2f:%.2f]\n", min, max);
}


// -------------------- gnuplot_one_function... ------------------------

void gnuplot_one_function_general (string title, string style, string label_x, string label_y, double *x, double *y, int n, bool jpg_file, string filename, bool square, bool ranges, double xmin, double xmax, double ymin, double ymax)
{
	int i;

	if (x == NULL || y == NULL || n < 1)
	{
		fprintf (stderr, "No data in gnuplot_one_function!\n");
		exit (EXIT_FAILURE);
	}

	gnuplot_init ();

	if (jpg_file)
	{
		fprintf (gnuplot_command, "set term jpeg\n\n\n");
		fprintf (gnuplot_command, "set out \"%s\"\n", filename.c_str ());
	}

	if (square)
		fprintf (gnuplot_command, "set size square\n");

	gnuplot_set_title (title);
	gnuplot_set_xlabel (label_x);
	gnuplot_set_ylabel (label_y);

	if (ranges)
	{
		gnuplot_set_xrange (xmin, xmax);
		gnuplot_set_yrange (ymin, ymax);
	}

	gnuplot_open_tmp_file ();
	for (i = 0; i < n; i++)
	{
		fprintf (gnuplot_tmp_file, "%.18e %.18e\n", x[i], y[i]);
	}
	gnuplot_close_tmp_file ();

	fprintf (gnuplot_command, "plot \"%s\" with %s\n", gnuplot_tmp_file_name.c_str (), style.c_str ());
	fflush (gnuplot_command);

	if (!jpg_file)
	{
		printf ("press ENTER to continue\n");
		while (getchar () != '\n')
		{
		}
		printf ("press ENTER to continue\n");
		while (getchar () != '\n')
		{
		}
	}
	gnuplot_close ();
	gnuplot_remove_tmp_file ();
}


void gnuplot_one_function (string title, string style, string label_x, string label_y, double *x, double *y, int n)
{
	gnuplot_one_function_general (title, style, label_x, label_y, x, y, n, false, "", false, false, 0.0, 0.0, 0.0, 0.0);
}


void gnuplot_one_function_jpg (string title, string style, string label_x, string label_y, double *x, double *y, int n, string filename)
{
	gnuplot_one_function_general (title, style, label_x, label_y, x, y, n, true, filename, false, false, 0.0, 0.0, 0.0, 0.0);
}


void gnuplot_one_function_square_jpg (string title, string style, string label_x, string label_y, double *x, double *y, int n, string filename)
{
	gnuplot_one_function_general (title, style, label_x, label_y, x, y, n, true, filename, true, false, 0.0, 0.0, 0.0, 0.0);
}


// ---------------- gnuplot_one_function_3d_... ------------------------

void gnuplot_one_function_3d_jpg (string title, string style, string label_x, string label_y, string label_z, double *x, double *y, void *z, int nx, int ny, int connected, double view1, double view2, string filename)
{
	int i, j;

	if (x == NULL || y == NULL || z == NULL || nx < 1 || ny < 1)
	{
		fprintf (stderr, "No data in gnuplot_one_function_3d!\n");
		exit (EXIT_FAILURE);
	}

	gnuplot_init ();

	if (filename != "")
	{
		fprintf (gnuplot_command, "set term jpeg\n\n\n");
		fprintf (gnuplot_command, "set out \"%s\"\n", filename.c_str ());
	}

	gnuplot_set_title (title);
	gnuplot_set_xlabel (label_x);
	gnuplot_set_ylabel (label_y);
	gnuplot_set_zlabel (label_z);

	fprintf (gnuplot_command, "set ticslevel 0\n");
	fprintf (gnuplot_command, "set view %f, %f\n", view1, view2);

	gnuplot_open_tmp_file ();
	for (i = 0; i < nx; i++)
	{
		for (j = 0; j < ny; j++)
			fprintf (gnuplot_tmp_file, "%.18e %.18e %.18e\n", x[i], y[j], ((double *) z)[j + ny * i]);
		fprintf (gnuplot_tmp_file, "\n");
		if (!connected)
			fprintf (gnuplot_tmp_file, "\n");
	}
	gnuplot_close_tmp_file ();

	fprintf (gnuplot_command, "splot \"%s\" using 1:2:3 with %s\n", gnuplot_tmp_file_name.c_str (), style.c_str ());
	fflush (gnuplot_command);

	if (filename == "")
	{
		printf ("press ENTER to continue\n");
		while (getchar () != '\n')
		{
		}
		printf ("press ENTER to continue\n");
		while (getchar () != '\n')
		{
		}
	}

	gnuplot_close ();
	gnuplot_remove_tmp_file ();
}


void gnuplot_one_function_3d (string title, string style, string label_x, string label_y, string label_z, double *x, double *y, double *z, int nx, int ny, int connected, double view1, double view2)
{
	gnuplot_one_function_3d_jpg (title, style, label_x, label_y, label_z, x, y, z, nx, ny, connected, view1, view2, "");
}


// ---------------- gnuplot_two_functions... ---------------------------


void gnuplot_two_functions_general (string title, string label_x, string label_y, double *x, double *y, int n, string dataname, string style, double *x2, double *y2, int n2, string dataname2, string style2, bool jpeg, string filename, bool square, bool ranges, double xmin, double xmax, double ymin, double ymax)
{
	int i;

	if (x == NULL || y == NULL || x2 == NULL || y2 == NULL || n < 1 || n2 < 1)
	{
		fprintf (stderr, "No data in gnuplot_two_functions!\n");
		exit (EXIT_FAILURE);
	}

	gnuplot_init ();

	if (jpeg)
	{
		fprintf (gnuplot_command, "set term jpeg\n\n\n");
		fprintf (gnuplot_command, "set out \"%s\"\n", filename.c_str ());
	}

	if (square)
		fprintf (gnuplot_command, "set size square\n");

	gnuplot_set_title (title);
	gnuplot_set_xlabel (label_x);
	gnuplot_set_ylabel (label_y);

	if (ranges)
	{
		gnuplot_set_xrange (xmin, xmax);
		gnuplot_set_yrange (ymin, ymax);
	}

	gnuplot_open_tmp_file ();
	for (i = 0; i < n; i++)
		fprintf (gnuplot_tmp_file, "%.18e %.18e\n", x[i], y[i]);
	fprintf (gnuplot_tmp_file, "\n\n");
	for (i = 0; i < n; i++)
		fprintf (gnuplot_tmp_file, "%.18e %.18e\n", x2[i], y2[i]);
	gnuplot_close_tmp_file ();

	fprintf (gnuplot_command, "plot \"%s\" index 0 title \"%s\" with %s, \"%s\" index 1 title \"%s\" with %s\n", gnuplot_tmp_file_name.c_str (), dataname.c_str (), style.c_str (), gnuplot_tmp_file_name.c_str (), dataname2.c_str (), style2.c_str ());
	fflush (gnuplot_command);

	if (!jpeg)
	{
		printf ("press ENTER to continue\n");
		while (getchar () != '\n')
		{
		}
		printf ("press ENTER to continue\n");
		while (getchar () != '\n')
		{
		}
	}
	gnuplot_close ();
	gnuplot_remove_tmp_file ();
}


void gnuplot_two_functions (string title, string style, string label_x, string label_y, double *x, double *y, int n, string dataname, double *x2, double *y2, int n2, string dataname2)
{
	gnuplot_two_functions_general (title, label_x, label_y, x, y, n, dataname, style, x2, y2, n2, dataname2, style, false, "", false, false, 0.0, 0.0, 0.0, 0.0);
}


void gnuplot_two_functions_jpg (string title, string style, string label_x, string label_y, double *x, double *y, int n, string dataname, double *x2, double *y2, int n2, string dataname2, string filename)
{
	gnuplot_two_functions_general (title, label_x, label_y, x, y, n, dataname, style, x2, y2, n2, dataname2, style, true, filename, false, false, 0.0, 0.0, 0.0, 0.0);
}


void gnuplot_two_functions_square_jpg (string title, string style, string label_x, string label_y, double *x, double *y, int n, string dataname, double *x2, double *y2, int n2, string dataname2, string filename)
{
	gnuplot_two_functions_general (title, label_x, label_y, x, y, n, dataname, style, x2, y2, n2, dataname2, style, true, filename, true, false, 0.0, 0.0, 0.0, 0.0);
}


// ----------------- Just create data file -------------------------------

void gnuplot_two_arrays_in_datafile (double *x, double *y, int n, string filename)
{
    FILE *file = fopen (filename.c_str(), "w");
    int nr;
    
    for (nr=0; nr< n; nr++)
        fprintf (file, "%e %e\n", x[nr], y[nr]);
    
    fclose (file);
}
