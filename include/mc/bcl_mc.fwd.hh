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

#ifndef BCL_MC_FWD_HH_
#define BCL_MC_FWD_HH_

// include the dependency file for this header
#include "bcl_mc.depends.fwd.hh"

// This file contains forward declarations for the mc namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace mc
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class MoviePrinterChimera;
    class MoviePrinterInterface;
    class MoviePrinterPymol;
    class MoviePrinters;
    class MutateLoopAdd;
    class MutateLoopAddResize;
    class MutateLoopFragmentAdd;
    class MutateLoopFragmentReplace;
    class MutateLoopRemove;
    class MutateLoopReplace;
    class Mutates;
    class OptimizationCCD;
    class OptimizationDocking;
    class OptimizationMCM;
    class Stage;
    class TemperatureAccepted;
    class TemperatureDefault;
    class TemperatureExponential;
    class TemperatureInterface;
    class TemperatureLinear;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class Approximator;

    template< typename t_ResultType>
    class Metropolis;

    template< typename t_ArgumentType, typename t_ResultType>
    class PrintInterface;

    template< typename t_ArgumentType, typename t_ResultType>
    class Printer;

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterCombined;

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterDefault;

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterFile;

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterWithCriterion;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< util::ShPtr< MoviePrinterInterface>, MoviePrinters> MoviePrinter;

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_FWD_HH_
