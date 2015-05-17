/*
 * =====================================================================================
 *
 *       Filename:  diagnostics.h
 *
 *    Description:  diagnostics for saving and reading data
 *
 *        Version:  1.0
 *        Created:  05/16/2015 13:34:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Di Qi (qidi), qidi@cims.nyu.edu
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _DIAGNOSTICS_H
#define _DIAGNOSTICS_H

// define file names
extern char ind_name[256];
extern char ug_name[256];
extern char vg_name[256];
extern char vorg_name[256];
extern char vors_name[256];

int diag_initial(char f_status);

int diag_writeData();
int diag_readData();

#endif /* DIAGNOSTICS_H */
