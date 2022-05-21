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
#include "fold/bcl_fold_mutate_sse_bend_template.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_fold_template.h"
#include "assemble/bcl_assemble_fold_template_handler.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateSSEBendTemplate::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateSSEBendTemplate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateSSEBendTemplate::MutateSSEBendTemplate() :
      m_SSEGeometryCompare( new assemble::SSEGeometryWithinSizeTolerance()),
      m_Scheme( GetStaticClassName< MutateSSEBendTemplate>())
    {
    }

    //! @brief construct from SSE geometry comparison function
    //! @param SSE_GEOMETRY_COMPARE SSE geometry comparison function
    //! @param SCHEME Scheme to be used
    MutateSSEBendTemplate::MutateSSEBendTemplate
    (
      const math::BinaryFunctionInterface< assemble::SSE, assemble::SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE,
      const std::string &SCHEME
    ) :
      m_SSEGeometryCompare( SSE_GEOMETRY_COMPARE.Clone()),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateSSEBendTemplate
    MutateSSEBendTemplate *MutateSSEBendTemplate::Clone() const
    {
      return new MutateSSEBendTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateSSEBendTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes an SSE and bends amino acids and returns the new SSE
    //! @param THIS_SSE SSE of interest to bend
    //! @return math::MutateResult that has a new SSE bent according to the random geometry from a fold template
    math::MutateResult< assemble::SSE> MutateSSEBendTemplate::operator()( const assemble::SSE &THIS_SSE) const
    {
      // initialize static undefined SSE
      static const util::ShPtr< assemble::SSE> sp_undefined_sse;

      // create a single element vector from the sse
      const util::SiPtr< const assemble::SSE> sp_sse( THIS_SSE);
      const util::SiPtrVector< const assemble::SSE> sse_vector( sp_sse);

      // get a random, single-geometry fold template
      const assemble::FoldTemplate fold_template
      (
        assemble::FoldTemplateHandler::GetRandomSubTemplate( sse_vector, *m_SSEGeometryCompare)
      );

      // if no fold template was found
      if( fold_template.GetGeometries().IsEmpty())
      {
        // return undefined mutate result
        return math::MutateResult< assemble::SSE>( sp_undefined_sse, *this);
      }

      // fit the sse
      const assemble::Domain fit_sse_domain( fold_template.FitSSEs( sse_vector));

      // if the domain is empty
      const util::SiPtrVector< const assemble::SSE> domain_sse_vector( fit_sse_domain.GetSSEs());
      if( domain_sse_vector.IsEmpty())
      {
        // return undefined mutate result
        return math::MutateResult< assemble::SSE>( sp_undefined_sse, *this);
      }

      // get the new sse
      util::ShPtr< assemble::SSE> new_sse( domain_sse_vector.FirstElement()->Clone());

      // move the new sse to the position of the original sse
      math::TransformationMatrix3D transform( math::Inverse( new_sse->GetOrientation()));
      transform( THIS_SSE.GetOrientation());
      new_sse->Transform( transform);

      // end
      return math::MutateResult< assemble::SSE>( new_sse, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateSSEBendTemplate::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSEGeometryCompare, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateSSEBendTemplate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSEGeometryCompare, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace fold

} // namespace bcl
