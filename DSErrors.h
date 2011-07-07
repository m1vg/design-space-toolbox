/**
 * \file DSErrors.h 
 * \brief Header file with functions for error and exception handling.
 *
 * \details
 * This file specifies the design space standard for error handling.
 * Contained here are the necessary macros and functions to succesfully report
 * the errors throughout the design space library.
 *
 * Copyright (C) 2010 Jason Lomnitz.\n\n
 *
 * This file is part of the Design Space Toolbox V2 (C Library).
 *
 * The Design Space Toolbox V2 is free software: you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * The Design Space Toolbox V2 is distributed in the hope that it will be 
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with the Design Space Toolbox. If not, see 
 * <http://www.gnu.org/licenses/>.
 *
 * \author Jason Lomnitz.
 * \date 2011
 */
#include <stdio.h>
#include <stdlib.h>
#include "DSTypes.h"

#ifndef __DS_ERRORS__
#define __DS_ERRORS__


/*********  Error Messages  ***********/

#define M_DS_NOFILE    "File not found"                //!< Message for no file found.
#define M_DS_NULL      "NULL pointer"                  //!< Message for NULL pointer.
#define M_DS_NOFORMAT  "Format not known"              //!< Message for unknown format.
#define M_DS_WRONG     "Inconsistent data"             //!< Message for inconsistent data being used.
#define M_DS_EXISTS    "Data already exists"           //!< Message for data aleady existing.
#define M_DS_NOTHREAD  "Thread not created"            //!< Message for no thread created.
#define M_DS_MALLOC    "Memory alloc failed"           //!< Message for failure to allocate data.
#define M_DS_NOT_IMPL  "Functionality not implemented" //!< Message for a feature not yet implemented.


/*********** Error Actions **************/

#define A_DS_NOERROR     0                          //!< Value for no error.
#define A_DS_WARN       -1                          //!< Value for a warning
#define A_DS_ERROR      -2                          //!< Value for an error.
#define A_DS_FATAL      -3                          //!< Value for a fatal error, kills program.
#define A_DS_KILLNOW    A_DS_FATAL                  //!< DEPRECATED: 


#define DSError(M_DS_Message, A_DS_Action) DSErrorFunction(M_DS_Message, A_DS_Action, __FILE__, __LINE__, __func__)

/**
 * \brief Pointer to a function determining how error messages are handled.
 *
 * This pointer to a function tells the error handling system which function to
 * call with the error messages.  If this pointer is NULL, then the system uses
 * printf, except that it prints to DSErrorFile instead of stdout.  This pointer
 * is intended to be used to override default behavior.  If used with MATLAB,
 * the pointer should be to mexPrintf.  If used in a Cocoa app, a function that
 * uses the notification system may be used.
 *
 * \see DSErrorFile
 * \see DSErrorFunction
 */
int (*DSPrintFunction)(const char *restrict, ...);

/**
 * \brief FILE pointer used for default DSPrintFunction.
 *
 * This pointer to a FILE tells the error handling system which FILE to
 * print the error messages to.  If this pointer is NULL, then the system uses
 * the stderr file.  This variable is only used internally with the default 
 * behavior of DSErrorFunction, however, it is intended to be used with 
 * custom functions.
 *
 * \see DSPrintFunction
 * \see DSErrorFunction
 */
FILE *DSErrorFile;

#ifdef __cplusplus
__BEGIN_DECLS
#endif


extern void DSErrorFunction(const char * M_DS_Message, char A_DS_ACTION, const char *FILEN, int LINE, const char * FUNC);



#ifdef __cplusplus
__END_DECLS    
#endif

#endif/* __DS_ERRORS__ */


/******************************* Documentation *******************************/
/**
 *\defgroup M_DS_Messages Messages for DS Errors.
 *
 * Defined here are the pre-defined messages used to report the
 * appropriate errors.  These are used with the different actions in
 * the macro DS_ERROR.  Other messages can be reported by literally
 * writting them in instead of these messages in the DS_ERROR macro.
 *
 *\see A_DS_Actions
 *\see DS_ERROR
 */
/*@{*/
/**
 * \def M_DS_NOFILE
 * \def M_DS_NULL
 * \def M_DS_NOFORMAT
 * \def M_DS_NOPARSE
 * \def M_DS_NOGMA
 * \def M_DS_NOSSYS
 */
/*@}*/

/**
 *\defgroup A_DS_Actions Actions for DS Errors.
 *
 * Defined here are the appropriate reactions to a specific error, an
 * error can have different actions depending on the sensitivity of
 * the region involved.
 *\see M_DS_Messages
 *\see DS_ERROR
 */
/*@{*/
/**
 * \def A_DS_NOERROR
 * \def A_DS_WARN
 * \def A_DS_ERROR
 * \def A_DS_KILLNOW
 */
/*@}*/

/**
 * \def DS_ERROR
 * \brief Error reporting macro.
 * \details
 *
 * Definition of the error reporting macro used within DS C
 * programs, this is a define which takes a message, some standard
 * messages are already included (see M_DS_Messages) and an action
 * (see A_DS_Actions) and reports to stderr.
 *
 * \see M_DS_Messages
 * \see A_DS_Actions
 */
