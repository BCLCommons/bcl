// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_H_
#define BCL_H_

// include the namespace forward header
#include "bcl.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <string>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
//! @mainpage BCL Documentation
//!
//! @image html BCL_m.jpg
//!
//! @section intro Introduction
//!
//! The BioChemistryLibrary (BCL) is a C++ library that includes classes to model
//! rather small molecules as organic substances as well as larger molecules as
//! proteins, dna, and rna. It contains forcefields, optimization algorithms,
//! as well as different prediction approaches as neural networks and support vector
//! machiens to model their structure, interactions, and properties.
//!
//! @author Meiler Group(
//! <A HREF="mailto: jens@jens-meiler.de"> jens@jens-meiler.de </A>
//! )\n
//! <A HREF="www.meilerlab.org"> www.meilerlab.org </A>\n
//!
//! @version 2.3
//!
//! @date 2002 - 2008
//!
//! @section rules Programming Rules
//!
//! @subsection variable Variable Name Conventions:
//! - global variables of namespaces are indicated by g_ and start with
//! one capital letter: g_Variable
//! - member variables of classes are indicatede by m_ and start with
//! one capital letter: m_Variable
//! - static variables of classes are indicatede by s_ and start with
//! one capital letter: s_Variable
//! - variables passed as function arguments are named in capital letters: VARIABLE
//! - local variables are named in small letters: variable
//!
//! @subsection general General Guidlines:
//! - code MUST be commented
//! - code MUST be commented using Doxygen style
//! - in particular each variable (local, member, or global), each function, and each class have to have an entry
//! - within functions separate blocks logically and add a comment in front of each block
//! - do not 'use namespace std' instead indicate every element used from a non-bcl
//! namespace (e.g. std) with std::
//! - use const as much as possible to indicate const arguments, pointers, or member
//! functions which do not alter the objects member variables! (the keyword 'mutable'
//! infront of a member variable allows changes althought const is used)
//! - Use const also for variables, that will not be changed within a scope, this will help the compiler to speed up compilation and runtime!
//! - avoid passing pointers to functions, instead derive a wrapper class from the
//! bcl::math::Function class
//! - to get an undefined state of type T1 use util::Undefined\< T1\>::s_Undefined (e.g. double a = util::GetUndefined\< double\>();)
//! to check whether a is undefined use IsDefined( a);
//! util::GetUndefined\< double\>(), util::GetUndefined\< size_t\>(), util::GetUndefined< int>() are defined for alternative usage in non-template classes
//! - to pass messages to a user use the BCL_Message( util::MessageLevel, Description) defined in util/bcl_util_message.h
//! and pass the appropiate util::MessageLevel (e_Silent, e_Critical, e_Standard, e_Verbose, e_Debug) and a description
//! - If you use util::Message( util::Message::e_Debug, Description), note that it will also output the filename, functionname and linenumber
//! - If you want to change the outputstream, please pass your stream to util::Message::SetOutputStream( your stream) - the default is std::cout
//! - for errors use BCL_Assert( Condition, Description) defined in util/bcl_util_assert.h. the program will terminate if Condition == false.
//! Pass an appropiate description. The Output will have, filename, function name and linenumber as a hint for error search
//! - every class should be derived from the ObjectInterface. Please read further descriptions in its documentation
//!
//! @subsection file File Organization:
//! - all code for a template class object must be included in the header (h). Put larger functions that
//! would usually go into the source code (cpp) file to the end (so that they can easily be transferred
//! later). This is currently done since the "export" keyword is not supported on template
//! functions in MVC. In turn it is impossible to put template class functions into separate
//! source code (cpp) files: The linker will not recognize the function specified for a particular
//! template class (error 2019). Given this fact the linker also reports multiple definitions of
//! in non-template classes (error 2005) for template functions of included template classes
//! if header (h) and source code (cpp) were separated! (This is rather anoying and no good style =\>
//! needs to be resolved)
//! - include an "example_CLASS.cpp" for every class defined in bcl_CLASS.h that demonstrates usage of the
//! class
//!
//! @subsection format util::Formating Guidlines:
//! - at least one space ' ' comes AFTER opening or closing a bracket but not if it is followed by "." or "->" or ")"
//! - no space comes after a statement befor closing a brace! like this: \n
//!     @code
//!     vector( 5);
//!     array[ 4];
//!     function( ARGUMENT1, ARGUMENT2);  @endcode
//!   this would be wrong:
//!     @code
//!     vector(7);
//!     vector(7 );
//!     vector( 7 );
//!     array[4];
//!     array[4 ];
//!     array[ 4 ];
//!     function( ARGUMENT1 , ARGUMENT2 );  @endcode
//!
//! - no spaces are left befor or after the '->' operator or "." operator
//! - indents are two space, not one tab, four spaces or whatever you like ( look in bcl code and you will see)
//! - blocks begin with '{' in the next line, for everything ( loops, declarations ...)
//! - blocks with one statement may be look like this: \n
//!   @code for( ...) @endcode
//!   { one statement;}
//! - blocks with multiple statements must look like this:
//!   @code for( ...)
//!   {
//!     statement1;
//!     statement1;
//!   }  @endcode
//! - This is a general rule: after each statement must be a new line!
//! - Do not use more than one empty line between statements, declarations, definitons, comments! \n
//!   @code function1( ...)
//!   {
//!     ...
//!   }
//!
//!   function2( ...)
//!   {
//!     variable a;
//
//!     ...
//!   }
//!
//!   //some comment
//!
//!   function3( ...)
//!   {
//!     ...
//!   }  @endcode
//!   \n
//!
//! @subsection ref_pointer References and Pointer
//! - ARGUMENTS shall be references unless otherwise inteted. Arguments that are not altered by the function have to be const!\n
//! - Memberfunctions that do not alter the objects have to be const member functions
//! - pointer to objects which are only used to acces the object without changing it should be const pointer. The same is true for util::SiPtr: util::SiPtr\< const type\>
//! - Functions that do take optional Arguments to return information to the user should have pointer to that argument, which are initilized with NULL as default.
//! - if the return is not optional, it is supposed to be a reference
//! - using references or pointer initializiation the '*' or '\&' comes infront of the variable name, and not as a part of the type name: \n
//!   right: @code
//!          double *a;
//!          SSE< AABackBone> *sse;
//!          ProteinModel< AACaCb> &PROTEINMODEL  @endcode
//!   wrong: @code
//!          double* a;
//!          SSE< AABackBone>* sse;
//!          ProteinModel< AACaCb>& PROTEINMODEL  @endcode
//!   This does not only count for local variables but also for function declarations!
//!
//! @section version SVN
//! - <A HREF="http://svnbook.red-bean.com/"> SVN online documentation </A>\n
//! - the directory of the bcl is on the accre-cluster: https://mimir.accre.vanderbilt.edu/repos/bcl
//!
//! @subsection download Checkout the current Version
//! - 'svn checkout --username="VUNET-ID" https://mimir.accre.vanderbilt.edu/repos/bcl'
//! - You will get a directory bcl. The current version of inculde, example and the doxygen documentation you will find within bcl/trunk/
//! - It is recommended not to modify the trunk - copy the trunk in directory: bcl/branch and work with this version \n
//!   you can also submit this branches to the repos on accre and use it. Finally you copy the modified files into trunk and commit this version
//! - SVN checkout with eclipse :  Remove .eclipse directory in your home directory before trying to setup SVN. \n
//!   go to SVN repository exploring perspective and click right button and select new -\> SVN repository location. \n
//!   enter https://mimir.accre.vanderbilt.edu/repos/bcl and follow the instructions.
//! - Checkout the trunk as a new C++ managed make project into your workspace. All project files are in the repository:
//! - the .cdtbuild, the .project and the .cdtproject are also checked out. You have to use the name bcl for your project in eclipse, otherwise the default compilation configurations will not work.
//! - for Visual studio there is a solution deposited and a README that helps you to set up you windows for the bcl.
//!
//! @subsection update Update your version
//! - There is no need to checkout a new version every day. If you want to make modifications, than go to your directory \n
//!   (trunk or branches/yourbranch) and type: 'svn update' and all your files will be updated, if there were any changes.
//! - befor you commit your changes, make an update and check that the bcl is still copiling and the examples are still running properly.
//! - do an update on a regular basis. If you do not, people may change code which may effect yours. If you do not update frequently you
//!   will get a lot of these changes later and it is hard to incooperate them all. Also commit your working classes frequently so that
//!   other programmers see what you have done, can use your code and do not have to write code again. Then they are also able to change
//!   your code if there are major changes and you do not get in trouble after one month without commits.
//!
//! @subsection upload Commit modified version
//! - got bcl/trunk or bcl/branch then type: 'svn commit' - if there were no changes from other people the modified files were \n
//!   transfered and they get a new revision number. Otherwise you got a message, that your files are out of date. Than type: 'svn update' \n
//!   if svn can, it will insert the changes to your file. Than you can commit it. \n
//!   Otherwise, you have to insert your changes to the current revisions of the files.
//!
//! @section gprof Profiling your code with GNU gprof
//! - <A HREF="http://www.gnu.org/software/binutils/manual/gprof-2.9.1/html_chapter/gprof_toc.html"> GNU gprof online documentation </A>\n
//! - Use ONLY if you want to have a detailed list of times and number of function calls in your code! Using gprof significantly reduces the speed of your executable.
//! - You will get two blocks of information as result.
//! - The first lists each function ordered by total used time. It also lists number of calls and time usage per call.
//! - The second is ordered hirarchical. It starts with the call of main() and tells how much time is used by the function itself and how much time by which child-function.
//!   The resulting file might be very large and you might have to scroll down pretty long to find it (tip: search for "index").
//!
//! @subsection gprof_eclipse How to use gprof in Eclipse
//! - in your project, go to: Preferences \> C++ Build \> GCC C++ Compiler \> Debugging:\n
//! set: "Generate gprof information (-pg)"
//! - the "-pg" option has to be set also for the linker:
//! GCC C++ Linker:\n Command: "g++ -pg"
//! - build project
//! - run project with all needed options, e.g. example_debug/example_debug.exe
//!   (this creates an unreadable gmon.out file)
//! - run gprof with the name of the executable and pipe output into a file:\n
//!   gprof example_debug/example_debug.exe \> times_examples.txt
//!
//! @subsection gprof_else How to use gprof in other environments
//! - set "-pg" option for both compiling and linking (see <A HREF="http://www.gnu.org/software/binutils/manual/gprof-2.9.1/html_chapter/gprof_toc.html"> documentation </A>)
//! - proceed as described above
//!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
//! @file bcl.h
//! @brief namespace header for everything in the biochemistry library
//!
//! @author meilerj
//! @date January, 2005
//! @namespace bcl global namespace of the biochemistry library
//! @see @link example_bcl.cpp @endlink
//!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace bcl
{

  //! @brief identifier for the name space
  //! @return the name of the namespace
  BCL_API
  const std::string &GetNamespaceIdentifier();

  //! @brief get the copyright
  //! @return the copyright
  BCL_API
  const std::string &GetCopyright();

  //! system dependent variables address size
  const size_t g_SizeOfAddress( sizeof( long *)); //!< size of pointer (4 for 32bit systems and 8 for 64bit systems)
  //! size of system in bit
  const size_t g_SizeOfSystem( 8 * g_SizeOfAddress); //!< size of system (32bit systems or 64bit systems)

} // namespace bcl

#endif //BCL_H_
