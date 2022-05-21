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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_pick_sse_short_loops.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PickSSEShortLoops::s_Instance
    (
      util::Enumerated< find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, DomainInterface> >::AddInstance( new PickSSEShortLoops)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @brief MAX_LOOP_LENGTH maximum number of residues between SSEs to be classified as short loop ( 5 by defult)
    PickSSEShortLoops::PickSSEShortLoops( const size_t MAX_LOOP_LENGTH) :
      m_MaxShortLoopLength( MAX_LOOP_LENGTH)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new PickSSEShortLoops which is a copy of this
    PickSSEShortLoops *PickSSEShortLoops::Clone() const
    {
      return new PickSSEShortLoops( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PickSSEShortLoops::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &PickSSEShortLoops::GetAlias() const
    {
      static const std::string s_alias( "PickSSEShortLoops");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PickSSEShortLoops::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Optimization implementation for Monte Carlo Metropolis algorithms.");
      serializer.AddInitializer
      (
        "max_short_loop_length",
        "maximum length of a loop to be considered short",
        io::Serialization::GetAgent( &m_MaxShortLoopLength)
      );

      return serializer;
    }
  ////////////////
  // operations //
  ////////////////

    //! @brief Picks one of SSEs with the Short Loops to the given domain
    //! @param SSE_LIST is the SiPtrList which provides the pool of assemble::SSE to pick from
    //! @param SSE_DOMAIN is the domain to which SSEs from the SSE_LIST will be compared for short loops
    //! @return returns SiPtr to the assemble::SSE object which has a short loops to one of the SSEs in SSE_DOMAIN
    util::SiPtr< const SSE>
    PickSSEShortLoops::Pick( const util::SiPtrList< const SSE> &SSE_LIST, const DomainInterface &SSE_DOMAIN) const
    {
      // get eligible list of sses by comparing to the protein model for short loops
      util::SiPtrList< const SSE> eligible_sses
      (
        SSE_DOMAIN.GetSSEsWithShortLoops( SSE_LIST, m_MaxShortLoopLength)
      );

      // none are found, return empty SiPtr
      if( eligible_sses.IsEmpty())
      {
        return util::SiPtr< const SSE>();
      }
      // if only one found return that one
      else if( eligible_sses.GetSize() == 1)
      {
        return eligible_sses.FirstElement();
      }
      // else more than one is found pick randomly
      else
      {
        return PickSSERandom().Pick( eligible_sses);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
