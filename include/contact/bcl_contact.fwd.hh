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

#ifndef BCL_CONTACT_FWD_HH_
#define BCL_CONTACT_FWD_HH_

// include the dependency file for this header
#include "bcl_contact.depends.fwd.hh"

// This file contains forward declarations for the contact namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace contact
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class AACorrelationFromFile;
    class CalculateCorrelationsMI;
    class CalculateCorrelationsSM;
    class CorrelationMatrix;
    class CorrelationStorageFile;
    class Data;
    class Map;
    class Order;
    class PredictionMap;
    class Recovery;
    class SSEPredictionMap;
    class Statistics;
    class TypeData;
    class Types;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_Member>
    class CalculateCorrelationsInterface;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< TypeData, Types> Type;

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_FWD_HH_
