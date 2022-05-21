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
#include "math/bcl_math_const_function.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_trigonometric_transition.h"
#include "score/bcl_score_restraint_xlink.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RestraintXlink::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new RestraintXlink())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from default values
    RestraintXlink::RestraintXlink() :
      m_Scheme( GetDefaultScheme()),
      m_Restraints(),
      m_ROGCalculator(),
      m_ScoringFunction()
    {
    }

    //! @brief construct from a restraint set and scheme
    //! @param SP_RESTRAINTS shared pointer to the restraints obtained from a cross-linking experiment
    //! @param TRANSITION_LENGTH length of the transition region of the scoring function
    //! @param SCHEME scheme of this score
    RestraintXlink::RestraintXlink
    (
      const util::ShPtr< util::ShPtrVector< restraint::AtomDistance> > &SP_RESTRAINTS,
      double TRANSITION_LENGTH,
      const std::string &SCHEME
    ) :
      m_Scheme( SCHEME),
      m_Restraints( SP_RESTRAINTS),
      m_ROGCalculator(),
      m_TransitionLength( TRANSITION_LENGTH),
      m_ScoringFunction( CreateScoringFunction( m_TransitionLength))
    {
    }

    //! @brief returns a pointer to a new RestraintXlink
    //! @return pointer to a new RestraintXlink
    RestraintXlink *RestraintXlink::Clone() const
    {
      return new RestraintXlink( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &RestraintXlink::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the scheme of this score
    //! @return the scheme of this score
    const std::string &RestraintXlink::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief returns the default scheme of this score
    //! @return the default scheme of this score
    const std::string &RestraintXlink::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "xlink");
      return s_default_scheme;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &RestraintXlink::GetAlias() const
    {
      static const std::string s_name( "RestraintXlink");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RestraintXlink::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      ( "scores the agreement of a protein model with cross-linking data");
      serializer.AddInitializer
      (
        "transition length",
        "The transition length used to call the scoring function",
        io::Serialization::GetAgent( &m_TransitionLength)
      );
      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief scores the agreement of a given protein model with data obtained from a cross-linking experiment
    //! @param PROTEIN_MODEL protein model for which to compute the agreement
    //! @return agreement score normalized for each restraint with -1 being the best and 0 being the worst agreement
    double RestraintXlink::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // sum up score for all restraints
      double sum_score( 0.0);
      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator it( m_Restraints->Begin()),
        it_end( m_Restraints->End());
        it != it_end;
        ++it
      )
      {
        // get the relevant data for evaluating the agreement of the model with this cross-link
        const restraint::AtomDistance &restraint( **it);
        const double xl_length( restraint.GetDistance()->GetDistance());
        const util::SiPtr< const biol::Atom> sp_atom_first( restraint.GetData().First()->LocateAtom( PROTEIN_MODEL));
        const util::SiPtr< const biol::Atom> sp_atom_second( restraint.GetData().Second()->LocateAtom( PROTEIN_MODEL));

        // if the end points of the restraint are not in the model, skip it
        if( !sp_atom_first.IsDefined() || !sp_atom_second.IsDefined())
        {
          continue;
        }

        // compute the agreement of the protein model with the given cross-link
        const double score( ComputeAgreement( PROTEIN_MODEL, *sp_atom_first, *sp_atom_second, xl_length));
        sum_score += score;
      }

      return sum_score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from an input stream
    //! @param ISTREAM input stream to read members from
    //! @return the input stream
    std::istream &RestraintXlink::Read( std::istream &ISTREAM)
    {
      // read members from input stream
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_Restraints, ISTREAM);
      io::Serialize::Read( m_ROGCalculator, ISTREAM);
      io::Serialize::Read( m_TransitionLength, ISTREAM);
      io::Serialize::Read( m_ScoringFunction, ISTREAM);

      return ISTREAM;
    }

    //! @brief writes members into an output stream
    //! @param OSTREAM output stream to write members into
    //! @INDENT number of indentations to use
    //! @return the output stream
    std::ostream &RestraintXlink::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members into output stream
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_ROGCalculator, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_TransitionLength, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_ScoringFunction, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool RestraintXlink::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_ScoringFunction = CreateScoringFunction( m_TransitionLength);
      return true;
    }

    //! @brief computes the agreement of the given protein model with the cross-link between the given atoms
    //! @param PROTEIN_MODEL protein model to compute the agreement for
    //! @param ATOM_FIRST first end point of the cross-link
    //! @param ATOM_SECOND second end point of the cross-link
    //! @param XL_LENGTH length of the cross-linker
    //! @return agreement of the protein model with the given cross-link
    double RestraintXlink::ComputeAgreement
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const biol::Atom &ATOM_FIRST,
      const biol::Atom &ATOM_SECOND,
      const double XL_LENGTH
    ) const
    {
      // compute opening angle defined by the end points of the cross-link
      const linal::Vector3D center( PROTEIN_MODEL.GetCenter());
      const linal::Vector3D cos_first( linal::UnitVector( center, ATOM_FIRST.GetCoordinates()));
      const linal::Vector3D cos_second( linal::UnitVector( center, ATOM_SECOND.GetCoordinates()));
      const double proj_angle( std::abs( linal::ProjAngle( cos_first, cos_second)));

      // approximate the distance that needs to be spanned by the cross-link by computing the perimeter of the circle
      const double radius( std::sqrt( m_ROGCalculator.SquareRadiusOfGyration( PROTEIN_MODEL)));
      const double perimeter( radius * proj_angle);
      const double length( perimeter);

      // evaluate the difference between cross-linker length and the length that must be spanned
      const double score( ( *m_ScoringFunction)( XL_LENGTH - length));

      return score;
    }

    //! @brief creates the function that scores the agreement of a protein model with cross-linking data
    //! @param TRANSITION_LENGTH length of the transition region
    //! @return shared pointer to the scoring function
    util::ShPtr< math::PiecewiseFunction> RestraintXlink::CreateScoringFunction( double TRANSITION_LENGTH)
    {
      // define the pieces of the function
      const math::Range< double> lower_range( -100.0, -TRANSITION_LENGTH / 2.0);
      const math::Range< double> trans_range
      (
        math::RangeBorders::e_LeftOpen,
        -TRANSITION_LENGTH / 2.0,
        TRANSITION_LENGTH / 2.0,
        math::RangeBorders::e_RightOpen
      );
      const math::Range< double> upper_range( TRANSITION_LENGTH / 2.0, 100.0);

      // define the function pieces
      const math::ConstFunction< double, double> lower_func( 0.0);
      const math::TrigonometricTransition trans_func( -TRANSITION_LENGTH / 2.0, TRANSITION_LENGTH / 2.0, 0.0, -1.0);
      const math::ConstFunction< double, double> upper_func( -1.0);

      // define the resulting function
      storage::List< storage::Pair< math::Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > > > args;
      args.Append
      (
        storage::Pair< math::Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
        (
          lower_range, util::CloneToShPtr( lower_func)
        )
      );
      args.Append
      (
        storage::Pair< math::Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
        (
          trans_range, util::CloneToShPtr( trans_func)
        )
      );
      args.Append
      (
        storage::Pair< math::Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
        (
          upper_range, util::CloneToShPtr( upper_func)
        )
      );
      util::ShPtr< math::PiecewiseFunction> sp_scoring_function( new math::PiecewiseFunction( args));

      return sp_scoring_function;
    }

  } // namespace score
} // namespace bcl
