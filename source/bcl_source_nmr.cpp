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
#include "nmr/bcl_nmr.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace nmr
} // namespace bcl
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
#include "nmr/bcl_nmr_rdc_container.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RDCContainer::s_Instance
    (
      GetObjectInstances().AddInstance( new RDCContainer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RDCContainer::RDCContainer() :
      m_ExperimentalValues(),
      m_CalculatedValues()
    {
    }

    //! @brief construct from RDC values
    //! @param EXP_VALUES experimental values
    //! @param CALC_VALUES calculated values
    RDCContainer::RDCContainer
    (
      const storage::Vector< double> &EXP_VALUES,
      const storage::Vector< double> &CALC_VALUES
    ) :
      m_ExperimentalValues( EXP_VALUES),
      m_CalculatedValues( CALC_VALUES)
    {
      BCL_Assert
      (
        m_ExperimentalValues.GetSize() == m_CalculatedValues.GetSize(),
        "# of experimental RDS (" + util::Format()( m_ExperimentalValues.GetSize()) + ") != # of calculated RDCs(" +
          util::Format()( m_CalculatedValues.GetSize()) + ")"
      );
    }

    //! @brief Clone function
    //! @return pointer to new RDCContainer
    RDCContainer *RDCContainer::Clone() const
    {
      return new RDCContainer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RDCContainer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RDCContainer::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExperimentalValues, ISTREAM);
      io::Serialize::Read( m_CalculatedValues, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RDCContainer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExperimentalValues, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CalculatedValues, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace nmr
  
} // namespace bcl
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
#include "nmr/bcl_nmr_residual_dipolar_coupling_least_square_deviation.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "linal/bcl_linal_matrix_inversion_moore_penrose.h"
#include "linal/bcl_linal_matrix_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
  //////////
  // data //
  //////////

    //! initialize single instance of that class
    const util::SiPtr< const util::ObjectInterface> ResidualDipolarCouplingLeastSquareDeviation::s_Instance
    (
      GetObjectInstances().AddInstance( new ResidualDipolarCouplingLeastSquareDeviation())
    );

    // initialize const size_t "number_independent_parameters" with 5
    const size_t ResidualDipolarCouplingLeastSquareDeviation::s_NumberIndependentParameters( 5);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ResidualDipolarCouplingLeastSquareDeviation::ResidualDipolarCouplingLeastSquareDeviation()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ResidualDipolarCoupling
    ResidualDipolarCouplingLeastSquareDeviation *ResidualDipolarCouplingLeastSquareDeviation::Clone() const
    {
      return new ResidualDipolarCouplingLeastSquareDeviation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ResidualDipolarCouplingLeastSquareDeviation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking a list of RDC assignments and returning a ResidualDipolarCouplingContainer
    //! The ResidualDipolarCouplingContainer will contain theoretical RDCs which have been calculated by doing a
    //! least squares deviation fitting to the experimental RDCs and the ResidualDipolarCouplingContainer will
    //! also contain the experimental RDCs
    //! @param RESTRAINTS list of RDC assignments
    //! @return experimental and calculated RDCs
    RDCContainer ResidualDipolarCouplingLeastSquareDeviation::operator()
    (
      const restraint::RDCAssignment &RESTRAINTS
    ) const
    {
      // create const size_t "number_of_datapoints" and initialize with the number of RDCs in the dataset
      const size_t number_datapoints( RESTRAINTS.GetData().GetSize());

      // if the number of restraints is too low to define a tensor
      if( number_datapoints < s_NumberIndependentParameters)
      {
        // return an empty container
        return RDCContainer();
      }

      // create math::Matrix of doubles "projection_angle_cosine"
      // initialize with the dimensions given by "number_datapoints" and "s_NumberIndependentParameters"
      linal::Matrix< double> projection_angle_cosine( number_datapoints, s_NumberIndependentParameters);

      // create math::Vector of doubles "experimental_rdcs" to hold the values of the experimental RDCs
      // initialize with the dimension given by "number_datapoints"
      linal::Vector< double> experimental_rdcs( number_datapoints);

      // create size_t "current_row" to access rows of "projection_angle_cosine" as the assignments
      // in "rdc_assignments" are iterated over in the for loop below
      // initialize with 0
      size_t current_row( 0);

      // iterate over the assignments to fill "projection_angle_cosine"
      for
      (
        storage::List< storage::Triplet< linal::Vector3D, linal::Vector3D, double> >::const_iterator
          itr( RESTRAINTS.GetData().Begin()), itr_end( RESTRAINTS.GetData().End());
        itr != itr_end;
        ++itr, ++current_row
      )
      {
        BCL_Assert( current_row < number_datapoints, "current row " + util::Format()( current_row));

        const double experimental_rdc_value( itr->Third());

        // assert these rather than continue because the matrix dimensions have already been set
        // it is up to another class to make sure at this point all the assignments are meaningful
        BCL_Assert
        (
          !math::EqualWithinAbsoluteTolerance( double( 0.0), experimental_rdc_value),
          "experimental_rdc_value is " + util::Format()( experimental_rdc_value)
        );

        // set the value of the current row of "experimental_rdcs" to the experimental RDC of "current_row"
        experimental_rdcs( current_row) = experimental_rdc_value;

        // create const linal::Vector3D "inter_nuclear_bond" to represent the bond between the two nuclei
        const linal::Vector3D inter_nuclear_bond( itr->First() - itr->Second());

        // create const double "cos_proj_ang_x_axis"
        const double cos_proj_ang_x_axis( std::cos( linal::ProjAngle( linal::Vector3D( 1, 0, 0), inter_nuclear_bond)));

        // create const double "cos_proj_ang_y_axis"
        const double cos_proj_ang_y_axis( std::cos( linal::ProjAngle( linal::Vector3D( 0, 1, 0), inter_nuclear_bond)));

        // create const double "cos_proj_ang_z_axis"
        const double cos_proj_ang_z_axis( std::cos( linal::ProjAngle( linal::Vector3D( 0, 0, 1), inter_nuclear_bond)));

        // create const double "sqr_cos_proj_ang_x_axis" and initialize with the square of "cos_proj_ang_x_axis"
        const double sqr_cos_proj_ang_x_axis( math::Sqr( cos_proj_ang_x_axis));

        const double double_cos_proj_ang_x_axis( double( 2) * cos_proj_ang_x_axis);

        // fill first column in current row of "projection_angle_cosine"
        projection_angle_cosine( current_row, 0) = math::Sqr( cos_proj_ang_y_axis) - sqr_cos_proj_ang_x_axis;

        // fill in second column in current row of "projection_angle_cosine"
        projection_angle_cosine( current_row, 1) = math::Sqr( cos_proj_ang_z_axis) - sqr_cos_proj_ang_x_axis;

        // fill in third column in current row of "projection_angle_cosine"
        projection_angle_cosine( current_row, 2) = double_cos_proj_ang_x_axis * cos_proj_ang_y_axis;

        // fill in fourth column in current row of "projection_angle_cosine"
        projection_angle_cosine( current_row, 3) = double_cos_proj_ang_x_axis * cos_proj_ang_z_axis;

        // fill in fifth column in current row of "projection_angle_cosine"
        projection_angle_cosine( current_row, 4) = double( 2) * cos_proj_ang_y_axis * cos_proj_ang_z_axis;
      }

      // invert "projection_angle_cosine"
      linal::MatrixInversionMoorePenrose< double> moore_penrose_inverter( projection_angle_cosine);

      // create math::Vector of doubles "independent_elements" and initialize with the solution to the
      // system of linear equations given by "projection_angle_cosine" and "experimental_rdcs"
      const linal::Vector< double> independent_elements( moore_penrose_inverter.Solve( experimental_rdcs));

      // create math Vector "theoretical_rdcs"
      // initialize with the RDCs as calculated from "not_inverted_projection_angle_cosine" and "independent_elements"
      const linal::Vector< double> theoretical_rdcs( projection_angle_cosine * independent_elements);

      // create pointer to double that is const "theoretical_rdcs_pointer"
      // initialize with the beginning of "theoretical_rdcs"
      const double *theoretical_rdcs_pointer( theoretical_rdcs.Begin());

      // create const pointer to double that is const "theoretical_rdcs_pointer_end"
      // initialize with the end of "theoretical_rdcs"
      const double *const theoretical_rdcs_pointer_end( theoretical_rdcs.End());

      // initialize rdc container
      storage::Vector< double> exp_rdcs;
      storage::Vector< double> calc_rdcs;

      // iterate over the restraints and store the calculated value
      for
      (
        storage::List< storage::Triplet< linal::Vector3D, linal::Vector3D, double> >::const_iterator
          rdc_itr( RESTRAINTS.GetData().Begin()), rdc_itr_end( RESTRAINTS.GetData().End());
        rdc_itr != rdc_itr_end && theoretical_rdcs_pointer != theoretical_rdcs_pointer_end;
        ++rdc_itr, ++theoretical_rdcs_pointer
      )
      {
        // pushback the values
        exp_rdcs.PushBack( rdc_itr->Third());
        calc_rdcs.PushBack( *theoretical_rdcs_pointer);
      }

      // return a ResidualDipolarCouplingContainer
      return RDCContainer( exp_rdcs, calc_rdcs);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ResidualDipolarCouplingLeastSquareDeviation::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ResidualDipolarCouplingLeastSquareDeviation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace nmr
} // namespace bcl
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
#include "nmr/bcl_nmr_rosetta_noe_handler.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"
#include "score/bcl_score_restraint_nmr_distance_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
  //////////
  // data //
  //////////

    //! defines the amount to be added to the NOE restraint for defining the upper limit
    const double RosettaNOEHandler::s_UpperLimit( 0.5);

    //! the lower limit for NOEs
    const double RosettaNOEHandler::s_LowerLimit( 1.8);

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RosettaNOEHandler::s_Instance
    (
      util::Enumerated< restraint::HandlerBase< util::ShPtrVector< restraint::AtomDistance> > >::AddInstance
      (
        new RosettaNOEHandler()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from the sequence distance and offset
    //! @param SEQUENCE_DISTANCE size_t which is the smallest distance in sequence two residues can be if a restraint
    //!        is going to be stored
    //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
    //! @param PREFIX prefix (path) for KB potential files
    //! @param SCORE_WEIGHT score weight to use
    RosettaNOEHandler::RosettaNOEHandler
    (
      const size_t SEQUENCE_DISTANCE,
      const size_t SEQUENCE_OFFSET,
      const std::string PREFIX,
      const double &SCORE_WEIGHT
    ) :
      m_SequenceDistance( SEQUENCE_DISTANCE),
      m_SequenceOffset( SEQUENCE_OFFSET),
      m_Prefix( PREFIX),
      m_Weight( SCORE_WEIGHT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RosettaNOEHandler
    RosettaNOEHandler *RosettaNOEHandler::Clone() const
    {
      return new RosettaNOEHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RosettaNOEHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &RosettaNOEHandler::GetAlias() const
    {
      static const std::string s_name( "Rosetta");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads restraints from a stream
    //! @param ISTREAM the stream the restraints will be read from
    //! @return the read in restraints
    util::ShPtrVector< restraint::AtomDistance> RosettaNOEHandler::ReadRestraints( std::istream &ISTREAM) const
    {
      // reset the restraints
      util::ShPtrVector< restraint::AtomDistance> restraints;

      // create a place to store the file
      storage::Vector< std::string> line_data( util::StringLineListFromIStream( ISTREAM));

      // iterate through the file's data
      for
      (
        storage::Vector< std::string>::const_iterator line_itr( line_data.Begin()), end_line_itr( line_data.End());
          line_itr != end_line_itr; ++line_itr
      )
      {
        // split the line
        storage::Vector< std::string> split_line( util::SplitString( *line_itr, "\t"));

        // break if the amount of data given is too small
        if( split_line.GetSize() < 10)
        {
          break;
        }

        // set chains to A
        const char chain_a( 'A');
        const char chain_b( 'A');

        // get the 1st and 2nd seq id and adjust according to the sequence offset
        const size_t seq_id_a( util::ConvertStringToNumericalValue< size_t>( split_line( 2)) - m_SequenceOffset);
        const size_t seq_id_b( util::ConvertStringToNumericalValue< size_t>( split_line( 4)) - m_SequenceOffset);

        // get the atom type strings
        const std::string atom_type_a( split_line( 1));
        const std::string atom_type_b( split_line( 3));

        // get the lower bound
        const double lower_bound( util::ConvertStringToNumericalValue< double>( split_line( 6)));

        // get the noe distance
        const double noe_file_value( util::ConvertStringToNumericalValue< double>( split_line( 7)));

        // get the actual noe value
        const double noe_value( std::max( s_LowerLimit, noe_file_value));

        // create an upper bound
        const double upper_bound( noe_value + s_UpperLimit);

        // create the distance restraint
        util::ShPtr< restraint::AtomDistance> sp_restraint
        (
          new restraint::AtomDistance
          (
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen( chain_a, seq_id_a, atom_type_a)
            ),
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen( chain_b, seq_id_b, atom_type_b)
            ),
            util::ShPtr< restraint::Distance>( new restraint::Distance( noe_value, upper_bound, lower_bound))
          )
        );

        // only store if if the restraint has defined values
        if( sp_restraint->IsDefined())
        {
          // only store if there is one chain and the two AAs are separated by enough amino acids
          if
          (
            math::Absolute( int( seq_id_a) - int( seq_id_b)) < int( m_SequenceDistance)
          )
          {
            // print a message
            BCL_MessageStd
            (
              "Sequence distance is smaller than cutoff, " +
              util::Format()( m_SequenceDistance) + " : " + sp_restraint->GetIdentification()
            );
          }
          // the sequence distance is large enough
          else
          {
            // add it to the vector
            restraints.PushBack( sp_restraint);
          }
        }
        else
        {
          // print a message
          BCL_MessageStd
          (
            "Rosetta file contains undefined restraint: " + sp_restraint->GetIdentification()
          );
        }
      }

      // end
      return restraints;
    }

    //! @brief writes restraints to a stream
    //! @param OSTREAM the stream the restraints will be written to
    //! @param RESTRAINT the restraint that will be written to the stream
    //! @return stream the restraints were written to
    std::ostream &RosettaNOEHandler::WriteRestraints
    (
      std::ostream &OSTREAM,
      const util::ShPtrVector< restraint::AtomDistance> &RESTRAINT
    ) const
    {
      // iterate through the vector of restraints
      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator itr( RESTRAINT.Begin()),
          end_itr( RESTRAINT.End());
        itr != end_itr; ++itr
      )
      {
        // dynamic cast the locators to LocatorCoordinatesHydrogen so that the string can be accessed
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_a( ( *itr)->GetData().First());
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_b( ( *itr)->GetData().Second());

        // assert that they are defined (otherwise these are not NOE restraints)
        BCL_Assert( sp_locator_a.IsDefined() && sp_locator_b.IsDefined(), "Restraints to be written are not NOEs");

        // print the type of ROSETTA restraint this is
        OSTREAM << "AtomPair\t";

        // if no prefix is given, use bounded potential
        if( m_Prefix == "")
        {
          // Give 1st Atom Type
          OSTREAM << sp_locator_a->GetAtomTypeString() << "\t";
          // Give 1st Atom SeqID
          OSTREAM << sp_locator_a->GetSeqID() << "\t";
          // Give 2nd Atom Type
          OSTREAM << sp_locator_b->GetAtomTypeString() << "\t";
          // Give 2nd Atom SeqID
          OSTREAM << sp_locator_b->GetSeqID() << "\t";
          // Required lines which tell how the restraint will be calculated (0) is the lower bound
          OSTREAM << "BOUNDED\t0\t";
          // and then the distance
          OSTREAM << ( *itr)->GetDistance()->GetDistance() << "\t";
          // tell that it's an NOE
          OSTREAM << util::Format()( m_Weight) << "\tNOE\n";
        }
        // use KB potentials
        else
        {
          // store the current atom types
          biol::AtomType atom_a( GetBaseAtomType( sp_locator_a->GetAtomTypeFromString()));
          biol::AtomType atom_b( GetBaseAtomType( sp_locator_b->GetAtomTypeFromString()));

          // get the total bonds to CB
          const size_t nr_bonds
          (
            score::RestraintNMRDistanceInterface::GetBondsFromCB( atom_a) +
            score::RestraintNMRDistanceInterface::GetBondsFromCB( atom_b)
          );

          // set the string
          std::string atom_string( "CB");
          if( atom_a == biol::GetAtomTypes().H || atom_b == biol::GetAtomTypes().H)
          {
            atom_string = "H";
          }
          if
          (
            atom_a == biol::GetAtomTypes().HA || atom_a == biol::GetAtomTypes().HA2 || atom_a == biol::GetAtomTypes().HA3 ||
            atom_b == biol::GetAtomTypes().HA || atom_b == biol::GetAtomTypes().HA2 || atom_b == biol::GetAtomTypes().HA3
          )
          {
            atom_string = "HA";
          }

          // handle HA2 and HA3
          std::string atom_a_name( atom_a.GetName());
          std::string atom_b_name( atom_b.GetName());
          if( atom_a == biol::GetAtomTypes().HA2)
          {
            atom_a_name = "1HA";
          }
          else if( atom_a == biol::GetAtomTypes().HA3)
          {
            atom_a_name = "2HA";
          }
          if( atom_b == biol::GetAtomTypes().HA2)
          {
            atom_b_name = "1HA";
          }
          else if( atom_b == biol::GetAtomTypes().HA3)
          {
            atom_b_name = "2HA";
          }

          // Give 1st Atom Type
          OSTREAM << atom_a_name << "\t";
          // Give 1st Atom SeqID
          OSTREAM << sp_locator_a->GetSeqID() << "\t";
          // Give 2nd Atom Type
          OSTREAM << atom_b_name << "\t";
          // Give 2nd Atom SeqID
          OSTREAM << sp_locator_b->GetSeqID() << "\t";
          // use bounded function if number of bonds is zero
          if( nr_bonds == 0)
          {
            // Required lines which tell how the restraint will be calculated (0) is the lower bound
            OSTREAM << "BOUNDED\t0\t";
            // and then the distance
            OSTREAM << ( *itr)->GetDistance()->GetDistance() << "\t";
            // tell that it's an NOE
            OSTREAM << util::Format()( m_Weight) << "\tNOE\n";
          }
          else
          {
            // output spline
            OSTREAM << "SPLINE\tdifference\t";
            // output filename
            OSTREAM << m_Prefix << atom_string << "_" << nr_bonds << ".potential" << '\t';
            // and then the distance
            OSTREAM << ( *itr)->GetDistance()->GetDistance() << "\t";
            OSTREAM << util::Format()( m_Weight) << "\t0.25\n";
          }
        }
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RosettaNOEHandler::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SequenceDistance, ISTREAM);
      io::Serialize::Read( m_SequenceOffset, ISTREAM);
      io::Serialize::Read( m_Prefix, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RosettaNOEHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SequenceDistance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SequenceOffset, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Prefix, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief gets the base atom type (i.e. CB)
    //! @param ATOM_TYPE atom type to be converted
    //! @return base atom type
    biol::AtomType RosettaNOEHandler::GetBaseAtomType( const biol::AtomType &ATOM_TYPE)
    {
      // initialize atom type to CB
      biol::AtomType atom_type( biol::GetAtomTypes().CB);

      // if type is HA or H
      if
      (
        ATOM_TYPE == biol::GetAtomTypes().HA ||
        ATOM_TYPE == biol::GetAtomTypes().HA2 ||
        ATOM_TYPE == biol::GetAtomTypes().HA3 ||
        ATOM_TYPE == biol::GetAtomTypes().H
      )
      {
        atom_type = ATOM_TYPE;
      }

      // end
      return atom_type;
    }

  } // namespace nmr
} // namespace bcl
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
#include "nmr/bcl_nmr_rosetta_rdc_handler.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RosettaRDCHandler::s_Instance
    (
      util::Enumerated< restraint::HandlerBase< restraint::RDC> >::AddInstance
      (
        new RosettaRDCHandler()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from the sequence offset
    //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
    RosettaRDCHandler::RosettaRDCHandler( const size_t SEQUENCE_OFFSET) :
      m_SequenceOffset( SEQUENCE_OFFSET),
      m_AdjustSigns( false),
      m_Normalize( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RosettaRDCHandler
    RosettaRDCHandler *RosettaRDCHandler::Clone() const
    {
      return new RosettaRDCHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RosettaRDCHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &RosettaRDCHandler::GetAlias() const
    {
      static const std::string s_name( "Rosetta");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads restraints from a stream
    //! @param ISTREAM the stream the restraints will be read from
    //! @return read in rdcs
    restraint::RDC RosettaRDCHandler::ReadRestraints( std::istream &ISTREAM) const
    {
      // reset the restraints
      restraint::RDC restraints;

      // store the file
      storage::Vector< std::string> line_data( util::StringLineListFromIStream( ISTREAM));

      // iterate through the line data
      for
      (
        storage::Vector< std::string>::const_iterator line_itr( line_data.Begin()), end_line_itr( line_data.End());
          line_itr != end_line_itr; ++line_itr
      )
      {
        // split each line by tabs
        storage::Vector< std::string> split_line( util::SplitString( *line_itr, "\t"));

        if( split_line.GetSize() < 5)
        {
          break;
        }

        // set chains to A
        const char chain_a( 'A');
        const char chain_b( 'A');

        // get the 1st and 2nd seq id and adjust according to the SequenceOffset
        const size_t seq_id_a( util::ConvertStringToNumericalValue< size_t>( split_line( 0)) - m_SequenceOffset);
        const size_t seq_id_b( util::ConvertStringToNumericalValue< size_t>( split_line( 2)) - m_SequenceOffset);

        // construct the locators
        const restraint::LocatorCoordinatesHydrogen locator_a( chain_a, seq_id_a, split_line( 1));
        const restraint::LocatorCoordinatesHydrogen locator_b( chain_b, seq_id_b, split_line( 3));

        // get the rdc value
        const double rdc_value( util::ConvertStringToNumericalValue< double>( split_line( 4)));

        // get the default bond length from the atom type
        const size_t bond_length( locator_a.GetAtomType()->GetBondLength( locator_b.GetAtomType()));

        // add the restraint to the restraint container
        restraints.PushBack( locator_a, locator_b, bond_length, rdc_value);
      } // end of looping through lines

      // if the normalize to NH flag is given
      if( m_Normalize)
      {
        // normalize the restraint
        restraints.NormalizetoNH();
      }
      // if the adjust sign flag is given
      else if( m_AdjustSigns)
      {
        // adjust the sign
        restraints.AdjustSigns();
      }

      // end
      return restraints;
    }

    //! @brief writes restraints to a stream
    //! @param OSTREAM the stream the restraints will be written to
    //! @param RESTRAINT the restraint that will be written to the stream
    //! @return stream the restraints were written to
    std::ostream &RosettaRDCHandler::WriteRestraints( std::ostream &OSTREAM, const restraint::RDC &RESTRAINT)
    {
      // iterate through the RDC container
      for
      (
        storage::Vector< storage::Triplet< restraint::DataPairwise, double, double> >::const_iterator
          rdc_itr( RESTRAINT.GetData().Begin()), end_itr( RESTRAINT.GetData().End());
        rdc_itr != end_itr; ++rdc_itr
      )
      {
        // dynamic cast the locators to LocatorCoordinatesHydrogen so that the string can be accessed
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_a( rdc_itr->First().First());
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_b( rdc_itr->First().Second());

        // assert that they are defined (otherwise these are not RDC restraints)
        BCL_Assert( sp_locator_a.IsDefined() && sp_locator_b.IsDefined(), "Restraints to be written are not RDCs");

        // seqID of first Atom
        OSTREAM << sp_locator_a->GetSeqID() << "\t";
        // type of first Atom
        OSTREAM << sp_locator_a->GetAtomTypeString() << "\t";
        // seqID of second Atom
        OSTREAM << sp_locator_b->GetSeqID() << "\t";
        // type of second Atom
        OSTREAM << sp_locator_b->GetAtomTypeString() << "\t";
        // RDC value
        OSTREAM << rdc_itr->Third() << "\n";
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RosettaRDCHandler::GetSerializer() const
    {
      io::Serializer serializer( restraint::HandlerInterface::GetSerializer());
      serializer.SetClassDescription
      (
        "Reads residual dipolar coupling restraints in Rosetta format; see "
        "http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/file_constraints.html for format specification"
      );
      serializer.AddInitializer
      (
        "offset",
        "offset of pdb ids relative to the restraints ids. This may be needed when the restraints were determined on a "
        "different sequence than is folded, defines the number to be subtracted from restraint ids to get to the "
        "sequence id ",
        io::Serialization::GetAgent( &m_SequenceOffset),
        "0"
      );
      serializer.AddInitializer
      (
        "normalize",
        "Normalize (sign and magnitude) of RDCs to NH values",
        io::Serialization::GetAgent( &m_Normalize),
        "False"
      );
      serializer.AddInitializer
      (
        "adjust signs",
        "adjust signs of normalized RDCs only. Set this if the values are scaled correctly, but non-N containing RDCs "
        "need to switch sign. Note that normalize=True will automatically adjust signs as well",
        io::Serialization::GetAgent( &m_AdjustSigns),
        "False"
      );
      return serializer;
    }

  } // namespace nmr
} // namespace bcl
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
#include "nmr/bcl_nmr_signal_1d.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

    //! @brief PatternType as string
    //! @param PATTERN_TYPE the PatternType
    //! @return the string for the PatternType
    const std::string &Signal1D::GetPatternTypeDescriptor( const PatternType &PATTERN_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "UnknownPattern",
        "Singlet",
        "Doublet",
        "Triplet",
        "Quadruplet",
        "Multiplet",
        GetStaticClassName< PatternType>()
      };

      return s_descriptors[ PATTERN_TYPE];
    }

    //! @brief PatternType as single character
    //! @details used in mdl Spectrum lines to describe the pattern
    //! @param PATTERN_TYPE the PatternType
    //! @return the char for the PatternType
    const char &Signal1D::GetPatternTypeSingleLetterCode( const PatternType &PATTERN_TYPE)
    {
      static const std::string s_single_letter_code( "USDTQMX");

      return s_single_letter_code[ PATTERN_TYPE];
    }

    //! @brief PatternType from single letter code
    //! @details convert char into a PatternType
    //! @param SINGLE_LETTER_CODE the single letter code as it may appear in an mdl Spectrum line
    //! @return PatternType for single letter code, e_UnknownPattern if not found
    Signal1D::PatternType Signal1D::PatternTypeFromSingleLetterCode( const char SINGLE_LETTER_CODE)
    {
      // range to search
      const char *beg( &GetPatternTypeSingleLetterCode( e_UnknownPattern));
      const char *end( beg + s_NumberPatternTypes);

      // find the single letter code
      const char *pattern( std::find( beg, end, SINGLE_LETTER_CODE));

      // if it was not found, return unknown
      if( pattern == end)
      {
        return e_UnknownPattern;
      }

      // return the PatternType
      return PatternType( pattern - beg);
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Signal1D::s_Instance
    (
      GetObjectInstances().AddInstance( new Signal1D())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Signal1D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Signal1D from std::istream
    std::istream &Signal1D::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ChemicalShift, ISTREAM);
      io::Serialize::Read( m_Integral,      ISTREAM);
      io::Serialize::Read( m_Pattern,       ISTREAM);
      io::Serialize::Read( m_Atom,          ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write Signal1D into std::ostream
    std::ostream &Signal1D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChemicalShift, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Integral,      OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Pattern,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Atom,          OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Signal from single entry in mdl
    //! @param SIGNAL_STRING string of the form shift;coupling_constcoupling;atom_index
    //! @param ATOMS the atoms, the atom index is relative to
    //! @return Signal1D for that string
    util::ShPtr< Signal1D> Signal1D::SignalFromString
    (
      const std::string &SIGNAL_STRING,
      iterate::Generic< const chemistry::AtomConformationalInterface> ATOMS
    )
    {
      // split the signal string
      const storage::Vector< std::string> signal_split( util::SplitString( SIGNAL_STRING, std::string( 1, s_InformationSeparator)));
      if( signal_split.GetSize() != 3)
      {
        BCL_MessageCrt
        (
          "each individual signal needs to contain: shift;coupling;atim_index, but it is: " + SIGNAL_STRING
        );
        return util::ShPtr< Signal1D>();
      }

      // shift and atom index
      const double current_shift( util::ConvertStringToNumericalValue< double>( signal_split( 0)));
      const size_t current_atom_index( util::ConvertStringToNumericalValue< size_t>( signal_split( 2)));

      // acquire the pattern
      // find the first character that is alpha
      const std::string::const_iterator first_pattern_itr( std::find_if( signal_split( 1).begin(), signal_split( 1).end(), std::ptr_fun( &isalpha)));
      const std::string coupling_const_string( signal_split( 1).begin(), first_pattern_itr);
      const std::string pattern_string( first_pattern_itr, signal_split( 1).end());

      const double coupling( util::IsNumerical( coupling_const_string) ? util::ConvertStringToNumericalValue< double>( coupling_const_string) : 0.0);
      const PatternType pattern( pattern_string.empty() ? e_UnknownPattern : PatternTypeFromSingleLetterCode( std::toupper( pattern_string[ pattern_string.length() - 1])));

      static const double s_default_integral( 1.0);

      ATOMS.GotoPosition( current_atom_index);
      // create signal 1D
      return util::ShPtr< Signal1D>
        (
          new Signal1D
          (
            current_shift,
            s_default_integral,
            storage::Vector< storage::Pair< PatternTypeEnum, double> >
            (
              1,
              storage::Pair< PatternTypeEnum, double>( pattern, coupling)
            ),
            *ATOMS
          )
        );
    }

  } // namespace nmr
} // namespace bcl

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
#include "nmr/bcl_nmr_signal.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Signal::s_Instance
    (
      GetObjectInstances().AddInstance( new Signal())
    );

  ////////////////
  // operations //
  ////////////////

    //! @brief check if given atom is involved in that signal
    //! @param ATOM the atom of interest
    //! @return true, if any Signal1D has the given atom
    bool Signal::ContainsAtom( const chemistry::AtomConformationalInterface &ATOM) const
    {
      const util::SiPtr< const chemistry::AtomConformationalInterface> si_ptr_atom( &ATOM);
      for
      (
        util::ShPtrVector< Signal1D>::const_iterator itr( m_Signals1D.Begin()), itr_end( m_Signals1D.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetAtomInvolvedInSignal() == si_ptr_atom)
        {
          return true;
        }
      }

      // no involved atoms agree
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Signal from std::istream
    //! @param ISTREAM input stream that contains Signal object
    std::istream &Signal::Read( std::istream &ISTREAM)
    {
      // read parameters
      io::Serialize::Read( m_Signals1D, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write Signal into std::ostream
    //! @param OSTREAM output stream into that the Signal object is written
    //! @param INDENT indentation
    std::ostream &Signal::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_Signals1D, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace nmr
} // namespace bcl

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
#include "nmr/bcl_nmr_spectrum.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

    //! @brief SpecType as string
    //! @param SPEC_TYPE the SpecType
    //! @return the string for the SpecType
    const std::string &Spectrum::GetSpecTypeDescriptor( const SpecType &SPEC_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "NO_SPECTYPE",
        "1H",
        "11B",
        "13C",
        "15N",
        "31P",
        "COSY",
        "DEPT",
        "HETCOR",
        "NOESY",
        GetStaticClassName< SpecType>()
      };

      return s_descriptors[ SPEC_TYPE];
    }

  //////////
  // data //
  //////////

    const std::string Spectrum::s_FieldStrengthPropertyDescriptor( "Field Strength [MHz]");
    const std::string Spectrum::s_TemperaturePropertyDescriptor( "Temperature [K]");
    const std::string Spectrum::s_SolventPropertyDescriptor( "Solvent");
    const std::string Spectrum::s_AssignmentMethodDescriptor( "Assignment Method");
    const std::string Spectrum::s_UnreportedDescriptor( "Unreported");
    const std::string Spectrum::s_SpectrumIdentifier( "Spectrum");

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Spectrum::s_Instance
    (
      GetObjectInstances().AddInstance( new Spectrum())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    Spectrum::Spectrum()
    {
    }

    //! @brief construct spectrum from its properties
    Spectrum::Spectrum( const SpecTypeEnum &TYPE, const util::ShPtrVector< Signal> &SIGNALS) :
      m_Type( TYPE),
      m_Signals( SIGNALS),
      m_Temperature( 298),
      m_Solvent(),
      m_FieldStrength( util::GetUndefined< double>())
    {
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief method for adding a spectrum to a small molecule
    //! @param MOLECULE the small molecule the spectrum is added to
    void Spectrum::AddSpectrum( chemistry::ConformationInterface &MOLECULE) const
    {
      std::string value;

      // for every signal
      for
      (
        util::ShPtrVector< Signal>::const_iterator
          itr_signal( m_Signals.Begin()), itr_signal_end( m_Signals.End());
        itr_signal != itr_signal_end;
        ++itr_signal
      )
      {
        const Signal1D &current_signal( *( *itr_signal)->GetSignals1D().FirstElement());

        // find the position of the atom that belongs to the given reference in signal1D object
        const size_t atom_index( MOLECULE.GetAtomIndex( *current_signal.GetAtomInvolvedInSignal()));

        const storage::Pair< Signal1D::PatternTypeEnum, double> &current_pattern( current_signal.GetPattern()( 0));

        // write spectrum    structure :  shift;couplingconstPattern;AtomIndex|
        value += util::Format().FFP( 1)( current_signal.GetChemicalShift()); // shift
        value += Signal1D::s_InformationSeparator;                   // separator
        value += util::Format().FFP( 1)( current_pattern.Second());  // coupling constant
        value += Signal1D::GetPatternTypeSingleLetterCode( current_pattern.First()); // coupling pattern
        value += Signal1D::s_InformationSeparator;                   // separator
        value += util::Format()( atom_index);                        // index of atom
        value += s_SignalSeparator;
      }

      // construct spectrum name from type and count
      const std::string spectrum_name_base( s_SpectrumIdentifier + " " + m_Type.GetString() + " ");

      // increase count until the spectrum does not exist as a misc property
      const storage::Map< size_t, util::ShPtr< Spectrum> > spectra( GenerateSpectra( MOLECULE));
      const size_t count( spectra.IsEmpty() ? 0 : spectra.ReverseBegin()->first + 1);

      // add the spectrum to the misc properties
      MOLECULE.StoreProperty( spectrum_name_base + util::Format()( count), value);

      // add properties

      // field strength
      storage::Pair< size_t, std::string> field_strength( count, s_UnreportedDescriptor);
      storage::Map< size_t, std::string> field_strength_props( GetPropertyMap( MOLECULE, s_FieldStrengthPropertyDescriptor));
      if( util::IsDefined( m_FieldStrength))
      {
        field_strength.Second() = util::Format()( m_FieldStrength);
      }
      BCL_Assert( field_strength_props.Insert( field_strength).second, "field strength already set: " + util::Format()( field_strength));
      MOLECULE.StoreProperty( s_FieldStrengthPropertyDescriptor, MapToLine( field_strength_props));

      // solvent
      storage::Map< size_t, std::string> solvent_props( GetPropertyMap( MOLECULE, s_SolventPropertyDescriptor));
      storage::Pair< size_t, std::string> solvent( count, s_UnreportedDescriptor);
      if( !m_Solvent.empty())
      {
        solvent.Second() = m_Solvent;
      }
      BCL_Assert( solvent_props.Insert( solvent).second, "solvent already set: " + util::Format()( solvent));
      MOLECULE.StoreProperty( s_SolventPropertyDescriptor, MapToLine( solvent_props));

      // assignment method
      storage::Map< size_t, std::string> assignment_method_props( GetPropertyMap( MOLECULE, s_AssignmentMethodDescriptor));
      storage::Pair< size_t, std::string> assignment( count, s_UnreportedDescriptor);
      if( !m_AssignmentMethod.empty())
      {
        assignment.Second() = m_AssignmentMethod;
      }
      BCL_Assert( assignment_method_props.Insert( assignment).second, "assignment already set: " + util::Format()( assignment));
      MOLECULE.StoreProperty( s_AssignmentMethodDescriptor, MapToLine( assignment_method_props));

      // temperature
      storage::Map< size_t, std::string> temperature_props( GetPropertyMap( MOLECULE, s_TemperaturePropertyDescriptor));
      storage::Pair< size_t, std::string> temperature( count, s_UnreportedDescriptor);
      if( util::IsDefined( m_Temperature))
      {
        temperature.Second() = util::Format()( m_Temperature);
      }
      BCL_Assert( temperature_props.Insert( temperature).second, "temperature already set: " + util::Format()( temperature));
      MOLECULE.StoreProperty( s_TemperaturePropertyDescriptor, MapToLine( temperature_props));
    }

    //! @brief get all signals involving given atom
    //! @param ATOM atom of interest
    //! @return SiPtrVector< const Signal> all signals involving atom of interest
    util::SiPtrVector< const Signal> Spectrum::DetermineInvolvedSignals( const chemistry::AtomConformationalInterface &ATOM) const
    {
      // all collected signals
      util::SiPtrVector< const Signal> signals;

      // iterate over all Signals in that Spectrum
      for
      (
        util::ShPtrVector< Signal>::const_iterator itr_signals( m_Signals.Begin()), itr_signals_end( m_Signals.End());
        itr_signals != itr_signals_end;
        ++itr_signals
      )
      {
        // check if current signal contains this atom
        if( ( *itr_signals)->ContainsAtom( ATOM))
        {
          signals.PushBack( util::SiPtr< const Signal>( **itr_signals));
        }
      }

      // return all signals for that atom
      return signals;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Spectrum from std::istream
    //! @param ISTREAM input stream that contains spectrum object
    std::istream &Spectrum::Read( std::istream &ISTREAM)
    {
      // read parameters
      io::Serialize::Read( m_Type            , ISTREAM);
      io::Serialize::Read( m_Signals         , ISTREAM);
      io::Serialize::Read( m_Temperature     , ISTREAM);
      io::Serialize::Read( m_Solvent         , ISTREAM);
      io::Serialize::Read( m_FieldStrength   , ISTREAM);
      io::Serialize::Read( m_AssignmentMethod, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write Spectrum into std::ostream
    //! @param OSTREAM output stream in which the spectrum object is written
    //! @param INDENT indentation
    std::ostream &Spectrum::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write parameters
      io::Serialize::Write( m_Type            , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Signals         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Temperature     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Solvent         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FieldStrength   , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AssignmentMethod, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief convert property line to property map
    //! @param CONFORMATION the conformation that should have the given property
    //! @param PROPERTY_LINE string of the property line "0:Unreported 1:CDCl3" create a map of size_t:string
    //! @return map that maps id to property
    storage::Map< size_t, std::string> Spectrum::GetPropertyMap
    (
      const chemistry::ConformationInterface &CONFORMATION,
      const std::string &PROPERTY_LINE
    )
    {
      const std::string &line( CONFORMATION.GetMDLProperty( PROPERTY_LINE));
      storage::Map< size_t, std::string> property_map;

      // split the string id separator
      const storage::Vector< std::string> split1
      (
        util::SplitString( line, std::string( 1, s_SpectraDescriptorSeparator))
      );

      if( split1.IsEmpty())
      {
        return property_map;
      }

      // get index of first element
      BCL_Assert
      (
        util::IsNumerical( split1.FirstElement()),
        "line needs to begin with an index, but is: " + line
      );
      size_t index( util::ConvertStringToNumericalValue< size_t>( split1.FirstElement()));
      for
      (
        storage::Vector< std::string>::const_iterator itr( split1.Begin() + 1), itr_last( split1.Last());
        itr <= itr_last;
        ++itr
      )
      {
        std::pair< size_t, std::string> key_value( index, "");

        // itr does not contain another index at the end
        if( itr == itr_last)
        {
          key_value.second = util::TrimString( *itr);
        }
        // itr contains the index for the next property at the end
        else
        {
          // find the property separator from the back
          const std::string::size_type split_pos( itr->rfind( s_PropertySeparator));
          BCL_Assert( split_pos != std::string::npos, "missing index for last property in line: " + line);

          // index is at the end
          const std::string next_index( itr->substr( split_pos, itr->length() - split_pos));
          BCL_Assert( util::IsNumerical( next_index), "non numerical index in line: " + line);
          index = util::ConvertStringToNumericalValue< size_t>( next_index);

          // actual property is before the index
          key_value.second = util::TrimString( itr->substr( 0, split_pos));
        }

        // check that it can be inserted without overlap of keys
        BCL_Assert( property_map.Insert( key_value).second, "properties with identical index: " + line);
      }

      // end
      return property_map;
    }

    //! @brief convert property map to property line "0:Unreported 1:CDCl3"
    //! @param PROPERTY_MAP that maps size_t to string
    //! @return string of the form "0:Unreported 1:CDCl3"
    std::string Spectrum::MapToLine( const storage::Map< size_t, std::string> &PROPERTY_MAP)
    {
      std::string property_line;

      // iterate over map
      for
      (
        storage::Map< size_t, std::string>::const_iterator itr( PROPERTY_MAP.Begin()), itr_end( PROPERTY_MAP.End());
        itr != itr_end;
        ++itr
      )
      {
        property_line += util::Format()( itr->first);
        property_line += s_SpectraDescriptorSeparator;
        property_line += itr->second;
        property_line += s_PropertySeparator;
      }

      // end
      return property_line;
    }

    //! @brief read signals from and MDL Spectrum line
    //! @param MOLECULE the small molecule the spectrum is for
    //! @param SPECTRUM_LINE the spectrum line with all the signals from mdl file misc properties in NMRSHIFTDB
    //! @return a SpHptrVector of Signals
    util::ShPtrVector< Signal> Spectrum::ReadSignals
    (
      const chemistry::ConformationInterface &MOLECULE,
      const std::string &SPECTRUM_LINE
    )
    {
      // shared pointer vector of signal1Ds
      util::ShPtrVector< Signal> vector_signal;

      // split spectrum line into signals
      const storage::Vector< std::string> signal_strings( util::SplitString( SPECTRUM_LINE, std::string( 1, s_SignalSeparator)));

      // iterate over all signal string
      for
      (
        storage::Vector< std::string>::const_iterator itr( signal_strings.Begin()), itr_end( signal_strings.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->length() < 2)
        {
          continue;
        }

        vector_signal.PushBack
        (
          util::ShPtr< Signal>
          (
            new Signal( util::ShPtrVector< Signal1D>( 1, Signal1D::SignalFromString( *itr, MOLECULE.GetAtomsIterator())))
          )
        );
      }

      // end
      return vector_signal;
    }

    //! @brief generate the spectra
    //! @param MOLECULE the small molecule associated with the spectra
    //! @return a Map of spectra, where the key is the spectra number
    storage::Map< size_t, util::ShPtr< Spectrum> > Spectrum::GenerateSpectra
    (
      const chemistry::ConformationInterface &MOLECULE
    )
    {
      // spectra properties
      const storage::Map< size_t, std::string> field_strength_props(    GetPropertyMap( MOLECULE, Spectrum::s_FieldStrengthPropertyDescriptor));
      const storage::Map< size_t, std::string> solvent_props(           GetPropertyMap( MOLECULE, Spectrum::s_SolventPropertyDescriptor));
      const storage::Map< size_t, std::string> temperature_props(       GetPropertyMap( MOLECULE, Spectrum::s_TemperaturePropertyDescriptor));
      const storage::Map< size_t, std::string> assignment_method_props( GetPropertyMap( MOLECULE, Spectrum::s_AssignmentMethodDescriptor));

      // initialize return value
      storage::Map< size_t, util::ShPtr< Spectrum> > spectra;

      // iterate over all properties and find Spectrum line
      for
      (
        storage::Map< std::string, std::string>::const_iterator
          prop_itr( MOLECULE.GetStoredProperties().Begin()), prop_itr_end( MOLECULE.GetStoredProperties().End());
        prop_itr != prop_itr_end;
        ++prop_itr
      )
      {
        const storage::Vector< std::string> descriptors( util::SplitString( prop_itr->first));
        if( descriptors( 0) != s_SpectrumIdentifier)
        {
          continue;
        }

        //BCL_MessageStd( MOLECULE.GetMiscellaneousProperty( "Spectrum 13C 0"));
        const size_t current_number( util::ConvertStringToNumericalValue< size_t>( descriptors( 2)));

        // check that spectra id is unique
        if( spectra.Find( current_number) != spectra.End())
        {
          BCL_MessageCrt( "spectra with duplicate id: " + descriptors( 2));
          continue;
        }

        // read out spectrum line in MDL block if available
        const util::ShPtrVector< Signal> signals
        (
          Spectrum::ReadSignals( MOLECULE, prop_itr->second)
        );

        // spectrum
        util::ShPtr< Spectrum> current_spectrum
        (
          new Spectrum( Spectrum::SpecTypeEnum( descriptors( 1)), signals)
        );

        // set temp and field strength, solvent and assignment method
        const storage::Map< size_t, std::string>::const_iterator temp_itr( temperature_props.Find(       current_number));
        const storage::Map< size_t, std::string>::const_iterator fiel_itr( field_strength_props.Find(    current_number));
        const storage::Map< size_t, std::string>::const_iterator solv_itr( solvent_props.Find(           current_number));
        const storage::Map< size_t, std::string>::const_iterator assi_itr( assignment_method_props.Find( current_number));

        // temperature
        if( temp_itr == temperature_props.End() || !util::IsNumerical( temp_itr->second))
        {
          // no temperature given
          current_spectrum->SetTemperature( 298.15);
        }
        else
        {
          current_spectrum->SetTemperature( util::ConvertStringToNumericalValue< double>( temp_itr->second));
        }

        // field strength
        if( fiel_itr == field_strength_props.End() || !util::IsNumerical( fiel_itr->second))
        {
          // no field strength given
          current_spectrum->SetFieldStrength( util::GetUndefined< double>());
        }
        else
        {
          current_spectrum->SetFieldStrength( util::ConvertStringToNumericalValue< double>( fiel_itr->second));
        }

        // solvent
        if( solv_itr == solvent_props.End() || solv_itr->second.empty())
        {
          current_spectrum->SetSolvent( Spectrum::s_UnreportedDescriptor);
        }
        else
        {
          current_spectrum->SetSolvent( solv_itr->second);
        }

        // assignment method
        if( assi_itr == assignment_method_props.End() || assi_itr->second.empty())
        {
          current_spectrum->SetAssignmentMethod( Spectrum::s_UnreportedDescriptor);
        }
        else
        {
          current_spectrum->SetAssignmentMethod( assi_itr->second);
        }

        // insert spectrum into map
        spectra[ current_number] = current_spectrum;
      }

      // end
      return spectra;
    }

  } // namespace nmr
} // namespace bcl
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
#include "nmr/bcl_nmr_star_noe_handler.h"

// includes from bcl - sorted alphabetically
#include "nmr/bcl_nmr_star_tags.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
  //////////
  // data //
  //////////

    //! defines the amount to be added to the NOE restraint for defining the upper limit
    const double StarNOEHandler::s_UpperLimit( 0.5);

    //! the lower limit for NOEs
    const double StarNOEHandler::s_LowerLimit( 1.8);

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> StarNOEHandler::s_Instance
    (
      util::Enumerated< restraint::HandlerBase< util::ShPtrVector< restraint::AtomDistance> > >::AddInstance
      (
        new StarNOEHandler()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from the sequence distance and offset
    //! @param SEQUENCE_DISTANCE size_t which is the smallest distance in sequence two residues can be if a restraint
    //!        is going to be stored
    //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
    StarNOEHandler::StarNOEHandler( const std::string &DEFAULT_EXTENSION, const size_t SEQUENCE_DISTANCE, const size_t SEQUENCE_OFFSET) :
      restraint::HandlerBase< util::ShPtrVector< restraint::AtomDistance> >( DEFAULT_EXTENSION),
      m_SequenceDistance( SEQUENCE_DISTANCE),
      m_SequenceOffset( SEQUENCE_OFFSET)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StarNOEHandler
    StarNOEHandler *StarNOEHandler::Clone() const
    {
      return new StarNOEHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarNOEHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &StarNOEHandler::GetAlias() const
    {
      static const std::string s_name( "Star");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads restraints from a stream
    //! @param ISTREAM the stream the restraints will be read from
    //! @return stream that the restraints were read from
    util::ShPtrVector< restraint::AtomDistance> StarNOEHandler::ReadRestraints( std::istream &ISTREAM) const
    {
      // reset the restraints
      util::ShPtrVector< restraint::AtomDistance> restraints;

      // read in the star file
      storage::Map
      <
        StarTagCategory,
        storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      > line_data( StarTags::ReadStarFile( ISTREAM));

      // first need to find the corresponding atom pairs from the distance constraint category

      // create a variable to hold the atom pair information
      storage::Map< size_t, storage::VectorND< 2, storage::Triplet< char, int, std::string> > > atom_pairs
      (
        GetAtomPairData( line_data)
      );

      // now read in the distance information that corresponds to stored atom pairs

      // create an iterator to the dist constraint info in the map
      storage::Map
      <
        StarTagCategory, storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      >::const_iterator map_itr_find( line_data.Find( GetStarTagCategories().e_DistConstraintValue));

      // make sure the input file has been read in and has dist constraint data
      BCL_Assert
      (
        map_itr_find != line_data.End(),
        "File containing distance restraints was not read in or does not contain the _Dist_constraint_value category"
      );

      // iterate through the dist constraint data
      for
      (
        storage::List< std::string>::const_iterator data_itr( map_itr_find->second.Second().Begin()),
          data_itr_end( map_itr_find->second.Second().End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // split the string by ' '
        storage::Vector< std::string> split_line( util::SplitString( *data_itr, " \t"));

        // initialize the distance, upper and lower bounds
        double distance;
        double upper_bound;
        double lower_bound;

        // if the distance is not in the correct place in the star file, define it by the upper bound
        if
        (
          split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceVal)) == "." ||
          split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceVal)) == "0.0"
        )
        {
          distance =
          util::ConvertStringToNumericalValue< double>
          (
            split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceUpperBoundVal))
          );
          upper_bound = distance + s_UpperLimit;

          // only give lower bound as standard if distance value is greater than 1.8
          if( distance > s_LowerLimit)
          {
            lower_bound = s_LowerLimit;
          }
          // else subtract the s_UpperLimit from the distance to get a lower bound
          else
          {
            lower_bound = distance - s_UpperLimit;
          }
        }

        // if the UpperBound value is missing, then calculate an upper and lower bound from the distance value
        else if
        (
          split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceUpperBoundVal)) == "."
        )
        {
          distance =
            util::ConvertStringToNumericalValue< double>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceVal))
            );
          upper_bound = distance + s_UpperLimit;
          // only give lower bound as standard if distance value is greater than 1.8
          if( distance > s_LowerLimit)
          {
            lower_bound = s_LowerLimit;
          }
          // else subtract the s_UpperLimit from the distance to get a lower bound
          else
          {
            lower_bound = distance - s_UpperLimit;
          }
        }
        else
        {
          distance =
            util::ConvertStringToNumericalValue< double>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceVal))
            );
          upper_bound =
            util::ConvertStringToNumericalValue< double>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceUpperBoundVal))
            );
          lower_bound =
            util::ConvertStringToNumericalValue< double>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueDistanceLowerBoundVal))
            );
        }

        // get the corresponding atom pair
        storage::VectorND< 2, storage::Triplet< char, int, std::string> > atom_pair
        (
          atom_pairs
          [
            util::ConvertStringToNumericalValue< size_t>
            (
              split_line( map_itr_find->second.First().Find( GetStarTags().e_DistValueConstraintID))
            )
          ]
        );

        const size_t seq_id_a( atom_pair.First().Second());
        const size_t seq_id_b( atom_pair.Second().Second());

        // create a new restraint
        util::ShPtr< restraint::AtomDistance> sp_restraint
        (
          new restraint::AtomDistance
          (
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen
              (
                atom_pair.First().First(),
                seq_id_a,
                atom_pair.First().Third()
              )
            ),
            util::ShPtr< assemble::LocatorAtomCoordinatesInterface>
            (
              new restraint::LocatorCoordinatesHydrogen
              (
                atom_pair.Second().First(),
                seq_id_b,
                atom_pair.Second().Third()
              )
            ),
            util::ShPtr< restraint::Distance>
            (
              new restraint::Distance
              (
                distance,
                upper_bound,
                lower_bound
              )
            )
          )
        );

        // if the restraint has defined values
        if( sp_restraint->IsDefined())
        {
          // if the sequence distance is too small and the atoms are on the same chain
          if
          (
            atom_pair.First().First() == atom_pair.Second().First() &&
            math::Absolute( atom_pair.First().Second() - atom_pair.Second().Second()) < int( m_SequenceDistance)
          )
          {
            // print a message
            BCL_MessageStd
            (
              "Sequence distance is smaller than cutoff, " +
              util::Format()( m_SequenceDistance) + " : " + sp_restraint->GetIdentification()
            );
          }
          // the sequence distance is large enough
          else
          {
            // add it to the vector
            restraints.PushBack( sp_restraint);
          }

        }
        // the restraint is undefined
        else
        {
          // print a message
          BCL_MessageStd
          (
            "Star file contains undefined restraint: " + sp_restraint->GetIdentification()
          );
        }
      }

      // end
      return restraints;
    }

    //! @brief writes restraints to a stream
    //! @param OSTREAM the stream the restraints will be written to
    //! @param RESTRAINT the restraint that will be written to the stream
    //! @return stream the restraints were written to
    std::ostream &StarNOEHandler::WriteRestraints
    (
      std::ostream &OSTREAM,
      const util::ShPtrVector< restraint::AtomDistance> &RESTRAINT
    )
    {
      // write file header
      const std::string dist_constraint( GetStarTagCategories().e_DistConstraint->GetDescription());
      OSTREAM << "loop_" << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistTreeNodeMemberConstraintID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistConstraintTreeNodeMemberID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistAuthAsymID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistAuthSeqID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_constraint << "." << GetStarTags().e_DistAuthAtomID->GetDescription() << '\n';

      size_t number_restraints( 0);

      // continue until the number of restraints has been satisfied or there are no more left to chose from
      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator itr( RESTRAINT.Begin()),
          end_itr( RESTRAINT.End());
        itr != end_itr; ++itr
      )
      {
        // dynamic cast the locators to LocatorCoordinatesHydrogen so that the string can be accessed
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_a( ( *itr)->GetData().First());
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_b( ( *itr)->GetData().Second());

        // assert that they are defined (otherwise these are not NOE restraints)
        BCL_Assert( sp_locator_a.IsDefined() && sp_locator_b.IsDefined(), "Restraints to be written are not NOEs");

        // make new line in output file
        OSTREAM << '\n';

        // write "Tree_node_member_constraint_ID"
        OSTREAM << number_restraints + 1 << ' ';
        // write "Constraint_tree_node_member_ID"
        OSTREAM << "1 ";
        // write Chain ID
        OSTREAM << sp_locator_a->GetChainID() << ' ';
        // write seq id
        OSTREAM << sp_locator_a->GetSeqID() << ' ';
        // write atom type
        OSTREAM << sp_locator_a->GetAtomTypeString();
        OSTREAM << '\n';
        // write "Tree_node_member_constraint_ID"
        OSTREAM << number_restraints + 1 << ' ';
        // write "Constraint_tree_node_member_ID"
        OSTREAM << "2 ";
        // write Chain ID
        OSTREAM << sp_locator_b->GetChainID() << ' ';
        // write seq id
        OSTREAM << sp_locator_b->GetSeqID() << ' ';
        // write atom type
        OSTREAM << sp_locator_b->GetAtomTypeString();

        ++number_restraints;
      }

      // write column headers
      const std::string dist_value( GetStarTagCategories().e_DistConstraintValue->GetDescription());
      OSTREAM << "\n\tstop_\n" << '\n';
      OSTREAM << "loop_" << '\n';
      OSTREAM << '\t' << dist_value << "." << GetStarTags().e_DistValueConstraintID->GetDescription() << '\n';
      OSTREAM << '\t' << dist_value << "." << GetStarTags().e_DistValueDistanceVal->GetDescription() << '\n';
      OSTREAM << '\t' << dist_value << "." << GetStarTags().e_DistValueDistanceLowerBoundVal->GetDescription() << '\n';
      OSTREAM << '\t' << dist_value << "." << GetStarTags().e_DistValueDistanceUpperBoundVal->GetDescription() << '\n';

      // loop through the vectors of AtomDistances to add to the bottom of the star file
      size_t counter( 1);

      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator itr( RESTRAINT.Begin()),
          itr_end( RESTRAINT.End());
        itr != itr_end; ++itr, ++counter
      )
      {
        OSTREAM << '\n';
        OSTREAM << counter;
        OSTREAM << ' ' << ( *itr)->GetDistance()->GetDistance();
        OSTREAM << ' ' << ( *itr)->GetDistance()->LowerBound();
        OSTREAM << ' ' << ( *itr)->GetDistance()->UpperBound();
      }

      OSTREAM << "\n\tstop_\nsave_";

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StarNOEHandler::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SequenceDistance, ISTREAM);
      io::Serialize::Read( m_SequenceOffset, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &StarNOEHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SequenceDistance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SequenceOffset, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief gets atom pair data used for determining distance restraints
    //! @param LINE_DATA the data needed to get NOE values
    //! @return vector of atom pair information (chain id, seq_id, atom type string)
    storage::Map
    <
      size_t, storage::VectorND< 2, storage::Triplet< char, int, std::string> >
    > StarNOEHandler::GetAtomPairData
    (
      storage::Map
      <
        StarTagCategory,
        storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      > &LINE_DATA
    ) const
    {
      // create an iterator to the dist constraint info in the map
      storage::Map
      <
        StarTagCategory, storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      >::const_iterator map_itr_find( LINE_DATA.Find( GetStarTagCategories().e_DistConstraint));

      // make sure the input file has been read in and has dist constraint data
      BCL_Assert
      (
        map_itr_find != LINE_DATA.End(),
        "File containing distance restraints was not read in or does not contain the _Dist_constraint category"
      );

      // create a variable to hold the atom pair information
      storage::Map< size_t, storage::VectorND< 2, storage::Triplet< char, int, std::string> > > atom_pairs;

      // create vector to hold chain ids
      storage::Map< std::string, char> chain_ids;

      // iterate through the dist constraint data
      for
      (
        storage::List< std::string>::const_iterator data_itr( map_itr_find->second.Second().Begin()),
          data_itr_end( map_itr_find->second.Second().End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // split the string by ' '
        storage::Vector< std::string> split_line_a( util::SplitString( *data_itr, " \t"));
        storage::Vector< std::string> split_line_b( util::SplitString( *( ++data_itr), " \t"));

        // get the current atom pair index and current atom
        const size_t current_index
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistTreeNodeMemberConstraintID))
          )
        );
        const size_t current_atom
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistConstraintTreeNodeMemberID))
          )
        );

        // create variable to hold the data
        storage::VectorND< 2, storage::Triplet< char, int, std::string> > atom_pair;

        // only read in the next line if it is a new atom (STAR format allows for multiple restraints from ambiguous
        // protons, such as 1   MET   HG2 vs 1   MET   HG3 that both may be in the restraint file, so only 1 needs to
        // be read in)
        while
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistTreeNodeMemberConstraintID))
          ) == current_index &&
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistConstraintTreeNodeMemberID))
          ) == current_atom
        )
        {
          ++data_itr;
          if( data_itr == data_itr_end)
          {
            break;
          }
          split_line_b = util::SplitString( *( data_itr), " \t");
        }

        // read in the chains
        const char chain_a
        (
          StarTags::GetChainID( split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistAuthAsymID)), chain_ids)
        );
        const char chain_b
        (
          StarTags::GetChainID( split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistAuthAsymID)), chain_ids)
        );

        // read in the seq ids
        const int seq_id_a
        (
          util::ConvertStringToNumericalValue< int>
          (
            split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistAuthSeqID))
          ) - int( m_SequenceOffset)
        );
        const int seq_id_b
        (
          util::ConvertStringToNumericalValue< int>
          (
            split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistAuthSeqID))
          ) - int( m_SequenceOffset)
        );

        BCL_Assert( util::IsDefined( m_SequenceDistance), "m_SequenceDistance is not defined");

        // read in the atom types as strings
        const std::string atom_a( split_line_a( map_itr_find->second.First().Find( GetStarTags().e_DistAuthAtomID)));
        const std::string atom_b( split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistAuthAtomID)));

        // add the atom pair to the vector
        storage::Triplet< char, int, std::string> atom_info_a( chain_a, seq_id_a, atom_a);
        storage::Triplet< char, int, std::string> atom_info_b( chain_b, seq_id_b, atom_b);
        atom_pair.First() = atom_info_a;
        atom_pair.Second() = atom_info_b;
        atom_pairs[ current_index] = atom_pair;

        // read through new lines if they are not part of a new atom pair
        while
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line_b( map_itr_find->second.First().Find( GetStarTags().e_DistTreeNodeMemberConstraintID))
          ) == current_index
        )
        {
          ++data_itr;
          if( data_itr == data_itr_end)
          {
            break;
          }
          split_line_b = util::SplitString( *( data_itr), " \t");
        }

        // found a new pair, so back up 1 so that it can be added in the next iteration
        --data_itr;
      }

      // end
      return atom_pairs;
    }

  } // namespace nmr
} // namespace bcl
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
#include "nmr/bcl_nmr_star_rdc_handler.h"

// includes from bcl - sorted alphabetically
#include "nmr/bcl_nmr_star_tags.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> StarRDCHandler::s_Instance
    (
      util::Enumerated< restraint::HandlerBase< restraint::RDC> >::AddInstance( new StarRDCHandler())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from the sequence offset
    //! @param SEQUENCE_OFFSET size_t which defines the number to be subtracted from the sequence id
    StarRDCHandler::StarRDCHandler
    (
      const std::string &DEFAULT_EXTENSION,
      const size_t SEQUENCE_OFFSET,
      const bool &ADJUST_SIGNS,
      const bool &NORMALIZE
    ) :
      restraint::HandlerBase< restraint::RDC>( DEFAULT_EXTENSION),
      m_SequenceOffset( SEQUENCE_OFFSET),
      m_AdjustSigns( ADJUST_SIGNS),
      m_Normalize( NORMALIZE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StarRDCHandler
    StarRDCHandler *StarRDCHandler::Clone() const
    {
      return new StarRDCHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarRDCHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns short name for the class used over the command line
    //! @return short name for the class used over the command line
    const std::string &StarRDCHandler::GetAlias() const
    {
      static const std::string s_name( "Star");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reads restraints from a stream
    //! @param ISTREAM the stream the restraints will be read from
    //! @return stream that the restraints were read from
    restraint::RDC StarRDCHandler::ReadRestraints( std::istream &ISTREAM) const
    {
      // reset the restraints
      restraint::RDC restraints;

      // read in the file.
      storage::Map
      <
        StarTagCategory,
        storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      > line_data( StarTags::ReadStarFile( ISTREAM));

      // find the RDC data in the map
      storage::Map
      <
        StarTagCategory, storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      >::const_iterator map_itr_find( line_data.Find( GetStarTagCategories().e_RDCConstraint));

      // make sure the input file has been read in and has RDC data
      BCL_Assert( map_itr_find != line_data.End(), "File containing RDC restraints was not read in");

      // create vector to hold chain ids
      storage::Map< std::string, char> chain_ids;

      // iterate through the rdc data lines
      for
      (
        storage::List< std::string>::const_iterator data_itr( map_itr_find->second.Second().Begin()),
          data_itr_end( map_itr_find->second.Second().End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // split the string by ' '
        storage::Vector< std::string> split_line( util::SplitString( *data_itr, " \t"));

        // create chains
        char chain_a
        (
          StarTags::GetChainID( split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthAsymID1)), chain_ids)
        );
        char chain_b
        (
          StarTags::GetChainID( split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthAsymID2)), chain_ids)
        );

        // create seq ids
        const size_t seq_id_a
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthSeqID1))
          ) - m_SequenceOffset
        );
        const size_t seq_id_b
        (
          util::ConvertStringToNumericalValue< size_t>
          (
            split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthSeqID2))
          ) - m_SequenceOffset
        );

        // create rdc value
        const double rdc_value
        (
          util::ConvertStringToNumericalValue< double>
          (
            split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCVal))
          )
        );

        // create atom types
        const std::string atom_type_a_string
        (
          split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthAtomID1))
        );
        const std::string atom_type_b_string
        (
          split_line( map_itr_find->second.First().Find( GetStarTags().e_RDCAuthAtomID2))
        );

        // construct the locators
        const restraint::LocatorCoordinatesHydrogen locator_a( chain_a, seq_id_a, atom_type_a_string);
        const restraint::LocatorCoordinatesHydrogen locator_b( chain_b, seq_id_b, atom_type_b_string);

        // find the bond length entry
        const size_t bond_length_index( map_itr_find->second.First().Find( GetStarTags().e_RDCBondLength));
        double bond_length;

        // if it was found
        if( bond_length_index != map_itr_find->second.First().GetSize())
        {
          // set the bond length
          bond_length = util::ConvertStringToNumericalValue< double>( split_line( bond_length_index));
        }
        // if it was not found
        else
        {
          // get the default bond length from the atom type
          bond_length = locator_a.GetAtomType()->GetBondLength( locator_b.GetAtomType());
          if( !util::IsDefined( bond_length))
          {
            continue;
          }
        }

        // add the restraint
        restraints.PushBack( locator_a, locator_b, bond_length, rdc_value);
      }

      // if the normalize to NH flag is given
      if( m_Normalize)
      {
        // normalize the restraint
        restraints.NormalizetoNH();
      }
      // if the adjust sign flag is given
      else if( m_AdjustSigns)
      {
        // adjust the sign
        restraints.AdjustSigns();
      }

      // end
      return restraints;
    }

    //! @brief writes restraints to a stream
    //! @param OSTREAM the stream the restraints will be written to
    //! @param RESTRAINT the restraint that will be written to the stream
    //! @return stream the restraints were written to
    std::ostream &StarRDCHandler::WriteRestraints( std::ostream &OSTREAM, const restraint::RDC &RESTRAINT)
    {
      // write file header
      const std::string rdc_constraint( GetStarTagCategories().e_RDCConstraint->GetDescription());
      OSTREAM << "loop_" << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCID->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCVal->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCLowerBound->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCUpperBound->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthAsymID1->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthSeqID1->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthAtomID1->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthAsymID2->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthSeqID2->GetDescription() << '\n';
      OSTREAM << '\t' << rdc_constraint << "." << GetStarTags().e_RDCAuthAtomID2->GetDescription() << '\n';

      // initialize a counter for writing out the NOE number
      int counter( 0);

      // iterate over the data
      for
      (
        storage::Vector< storage::Triplet< restraint::DataPairwise, double, double> >::const_iterator
          rdc_itr( RESTRAINT.GetData().Begin()), end_itr( RESTRAINT.GetData().End());
        rdc_itr != end_itr; ++rdc_itr, ++counter
      )
      {
        // dynamic cast the locators to LocatorCoordinatesHydrogen so that the string can be accessed
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_a( rdc_itr->First().First());
        const util::ShPtr< restraint::LocatorCoordinatesHydrogen> sp_locator_b( rdc_itr->First().Second());

        // assert that they are defined (otherwise these are not RDC restraints)
        BCL_Assert( sp_locator_a.IsDefined() && sp_locator_b.IsDefined(), "Restraints to be written are not RDCs");

        // print out data to file
        OSTREAM << '\n' << counter + 1 << ' ' << util::Format().FFP( 1)( rdc_itr->Third()) << ' '; // counter and val
        OSTREAM << util::Format().FFP( 1)( rdc_itr->Third() - 1) << ' '; // lower bound
        OSTREAM << util::Format().FFP( 1)( rdc_itr->Third() + 1) << ' '; // upper bound
        OSTREAM << sp_locator_a->GetChainID() << ' '; // chain id 1
        OSTREAM << sp_locator_a->GetSeqID() << ' '; // seq id 1
        OSTREAM << sp_locator_a->GetAtomTypeString() << ' '; // atom type 1
        OSTREAM << sp_locator_b->GetChainID() << ' '; // chain id 2
        OSTREAM << sp_locator_b->GetSeqID() << ' '; // seq id 2
        OSTREAM << sp_locator_b->GetAtomTypeString(); // atom type 2
      }

      OSTREAM << "\n\tstop_\nsave_";

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer StarRDCHandler::GetSerializer() const
    {
      io::Serializer serializer( restraint::HandlerInterface::GetSerializer());
      serializer.SetClassDescription
      (
        "Reads residual dipolar coupling restraints in STAR format; see "
        "http://www.bmrb.wisc.edu/dictionary/3.1html/SuperGroupPage.html for format specification"
      );
      serializer.AddInitializer
      (
        "offset",
        "offset of pdb ids relative to the restraints ids. This may be needed when the restraints were determined on a "
        "different sequence than is folded, defines the number to be subtracted from restraint ids to get to the "
        "sequence id ",
        io::Serialization::GetAgent( &m_SequenceOffset),
        "0"
      );
      serializer.AddInitializer
      (
        "normalize",
        "Normalize (sign and magnitude) of RDCs to NH values",
        io::Serialization::GetAgent( &m_Normalize),
        "False"
      );
      serializer.AddInitializer
      (
        "adjust signs",
        "adjust signs of normalized RDCs only. Set this if the values are scaled correctly, but non-N containing RDCs "
        "need to switch sign. Note that normalize=True will automatically adjust signs as well",
        io::Serialization::GetAgent( &m_AdjustSigns),
        "False"
      );
      return serializer;
    }

  } // namespace nmr
} // namespace bcl
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
#include "nmr/bcl_nmr_star_tag_categories.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    StarTagCategories::StarTagCategories() :
      e_DistConstraint(      AddEnum( "DistConstraint",      StarTagCategoryData( StarTagCategoryData::e_DistanceConstraints, "_Dist_constraint"))),
      e_DistConstraintList(  AddEnum( "DistConstraintList",  StarTagCategoryData( StarTagCategoryData::e_DistanceConstraints, "_Dist_constraint_list"))),
      e_DistConstraintTree(  AddEnum( "DistConstraintTree",  StarTagCategoryData( StarTagCategoryData::e_DistanceConstraints, "_Dist_constraint_tree"))),
      e_DistConstraintValue( AddEnum( "DistConstraintValue", StarTagCategoryData( StarTagCategoryData::e_DistanceConstraints, "_Dist_constraint_value"))),
      e_RDCConstraint(       AddEnum( "RDCConstraint",       StarTagCategoryData( StarTagCategoryData::e_RDCConstraints,      "_RDC_constraint"))),
      e_RDCConstraintList(   AddEnum( "RDCConstraintList",   StarTagCategoryData( StarTagCategoryData::e_RDCConstraints,      "_RDC_constraint_list")))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarTagCategories::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a tag category from a string
    //! @param DESCRIPTION string description for the enum
    //! @return a tag category from a string
    StarTagCategory StarTagCategories::GetCategoryFromString( const std::string &DESCRIPTION)
    {
      // iterate through the tag categories
      for
      (
        StarTagCategories::const_iterator category_itr( GetStarTagCategories().Begin()),
          category_itr_end( GetStarTagCategories().End());
        category_itr != category_itr_end;
        ++category_itr
      )
      {
        // if the category has the same description
        if( ( *category_itr)->GetDescription() == DESCRIPTION)
        {
          // return it
          return *category_itr;
        }
      }

      // no matching description found
      return GetStarTagCategories().e_Undefined;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief get access to all star tag categories
    const StarTagCategories &GetStarTagCategories()
    {
      return StarTagCategories::GetEnums();
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace nmr

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< nmr::StarTagCategoryData, nmr::StarTagCategories>;

  } // namespace util
} // namespace bcl
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
#include "nmr/bcl_nmr_star_tag_category_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

    //! @brief SaveFrame as string
    //! @param FRAME the save frame
    //! @return the string for FRAME
    const std::string &StarTagCategoryData::GetSaveFrameName( const SaveFrame &FRAME)
    {
      static const std::string s_names[] =
      {
        "DistanceConstraint",
        "RDCConstraint",
        GetStaticClassName< SaveFrame>()
      };
      return s_names[ FRAME];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from StarSaveFrame
    //! @param SAVE_FRAME save frame to assign to this tag category
    //! @param DESCRIPTION string used to recognize this category in an NMR-STAR file
    StarTagCategoryData::StarTagCategoryData( const SaveFrame &SAVE_FRAME, const std::string &DESCRIPTION) :
      m_SaveFrame( SAVE_FRAME),
      m_Description( DESCRIPTION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StarTagCategoryData
    StarTagCategoryData *StarTagCategoryData::Clone() const
    {
      return new StarTagCategoryData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarTagCategoryData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StarTagCategoryData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SaveFrame, ISTREAM);
      io::Serialize::Read( m_Description, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &StarTagCategoryData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SaveFrame, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace nmr
} // namespace bcl
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
#include "nmr/bcl_nmr_star_tag_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from tag category and data type
    //! @param CATEGORY category for this tag
    //! @param DESCRIPTION description string for this tag
    //! @param DATA_TYPE data type for this tag
    StarTagData::StarTagData
    (
      const StarTagCategory &CATEGORY,
      const std::string &DESCRIPTION,
      const util::CPPDataTypes::Types DATA_TYPE
    ) :
      m_TagCategory( CATEGORY),
      m_Description( DESCRIPTION),
      m_DataType( DATA_TYPE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StarTagData
    StarTagData *StarTagData::Clone() const
    {
      return new StarTagData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarTagData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StarTagData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_TagCategory, ISTREAM);
      io::Serialize::Read( m_Description, ISTREAM);
      io::Serialize::Read( m_DataType, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &StarTagData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_TagCategory, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DataType, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace nmr
} // namespace bcl
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
#include "nmr/bcl_nmr_star_tags.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace nmr
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    StarTags::StarTags() :
      e_DistTreeConstraintID(              AddEnum( "DistTreeConstraintID",              StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Constraint_ID",                  util::CPPDataTypes::e_SizeT))),
      e_DistTreeNodeID(                    AddEnum( "DistTreeNodeID",                    StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Node_ID",                        util::CPPDataTypes::e_SizeT))),
      e_DistTreeDownNodeID(                AddEnum( "DistTreeDownNodeID",                StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Down_node_ID",                   util::CPPDataTypes::e_SizeT))),
      e_DistTreeRightNodeID(               AddEnum( "DistTreeRightNodeID",               StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Right_node_ID",                  util::CPPDataTypes::e_SizeT))),
      e_DistTreeLogicOperation(            AddEnum( "DistTreeLogicOperation",            StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Logic_operation",                util::CPPDataTypes::e_String))),
      e_DistTreeEntryID(                   AddEnum( "DistTreeEntryID",                   StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_DistTreeDistanceConstraintListID(  AddEnum( "DistTreeDistanceConstraintListID",  StarTagData( GetStarTagCategories().e_DistConstraintTree,  "Distance_constraint_list_ID",    util::CPPDataTypes::e_SizeT))),

      e_DistTreeNodeMemberConstraintID(    AddEnum( "DistTreeNodeMemberConstraintID",    StarTagData( GetStarTagCategories().e_DistConstraint,      "Tree_node_member_constraint_ID", util::CPPDataTypes::e_SizeT))),
      e_DistTreeNodeMemberNodeID(          AddEnum( "DistTreeNodeMemberNodeID",          StarTagData( GetStarTagCategories().e_DistConstraint,      "Tree_node_member_node_ID",       util::CPPDataTypes::e_SizeT))),
      e_DistConstraintTreeNodeMemberID(    AddEnum( "DistConstraintTreeNodeMemberID",    StarTagData( GetStarTagCategories().e_DistConstraint,      "Constraint_tree_node_member_ID", util::CPPDataTypes::e_SizeT))),
      e_DistAssemblyAtomID(                AddEnum( "DistAssemblyAtomID",                StarTagData( GetStarTagCategories().e_DistConstraint,      "Assemble_atom_ID",               util::CPPDataTypes::e_SizeT))),
      e_DistEntityAssemblyID(              AddEnum( "DistEntityAssemblyID",              StarTagData( GetStarTagCategories().e_DistConstraint,      "Entity_assembly_ID",             util::CPPDataTypes::e_SizeT))),
      e_DistEntityID(                      AddEnum( "DistEntityID",                      StarTagData( GetStarTagCategories().e_DistConstraint,      "Entity_ID",                      util::CPPDataTypes::e_SizeT))),
      e_DistCompIndexID(                   AddEnum( "DistCompIndexID",                   StarTagData( GetStarTagCategories().e_DistConstraint,      "Comp_index_ID",                  util::CPPDataTypes::e_SizeT))),
      e_DistSeqID(                         AddEnum( "DistSeqID",                         StarTagData( GetStarTagCategories().e_DistConstraint,      "Seq_ID",                         util::CPPDataTypes::e_SizeT))),
      e_DistCompID(                        AddEnum( "DistCompID",                        StarTagData( GetStarTagCategories().e_DistConstraint,      "Comp_ID",                        util::CPPDataTypes::e_String))),
      e_DistAtomID(                        AddEnum( "DistAtomID",                        StarTagData( GetStarTagCategories().e_DistConstraint,      "Atom_ID",                        util::CPPDataTypes::e_String))),
      e_DistResonanceID(                   AddEnum( "DistResonanceID",                   StarTagData( GetStarTagCategories().e_DistConstraint,      "Resonance_ID",                   util::CPPDataTypes::e_SizeT))),
      e_DistAuthAsymID(                    AddEnum( "DistAuthAsymID",                    StarTagData( GetStarTagCategories().e_DistConstraint,      "Auth_asym_ID",                   util::CPPDataTypes::e_String))),
      e_DistAuthSeqID(                     AddEnum( "DistAuthSeqID",                     StarTagData( GetStarTagCategories().e_DistConstraint,      "Auth_seq_ID",                    util::CPPDataTypes::e_String))),
      e_DistAuthCompID(                    AddEnum( "DistAuthCompID",                    StarTagData( GetStarTagCategories().e_DistConstraint,      "Auth_comp_ID",                   util::CPPDataTypes::e_String))),
      e_DistAuthAtomID(                    AddEnum( "DistAuthAtomID",                    StarTagData( GetStarTagCategories().e_DistConstraint,      "Auth_atom_ID",                   util::CPPDataTypes::e_String))),
      e_DistEntryID(                       AddEnum( "DistEntryID",                       StarTagData( GetStarTagCategories().e_DistConstraint,      "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_DistDistanceConstraintListID(      AddEnum( "DistDistanceConstraintListID",      StarTagData( GetStarTagCategories().e_DistConstraint,      "Distance_constraint_list_ID",    util::CPPDataTypes::e_SizeT))),

      e_DistListSfCatetory(                AddEnum( "DistListSfCatetory",                StarTagData( GetStarTagCategories().e_DistConstraintList,  "Sf_category",                    util::CPPDataTypes::e_String))),
      e_DistListEntryID(                   AddEnum( "DistListEntryID",                   StarTagData( GetStarTagCategories().e_DistConstraintList,  "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_DistListSfID(                      AddEnum( "DistListSfID",                      StarTagData( GetStarTagCategories().e_DistConstraintList,  "Sf_ID",                          util::CPPDataTypes::e_SizeT))),
      e_DistListID(                        AddEnum( "DistListID",                        StarTagData( GetStarTagCategories().e_DistConstraintList,  "ID",                             util::CPPDataTypes::e_SizeT))),
      e_DistListConstraintType(            AddEnum( "DistListConstraintType",            StarTagData( GetStarTagCategories().e_DistConstraintList,  "Constraint_type",                util::CPPDataTypes::e_String))),
      e_DistListConstraintFileID(          AddEnum( "DistListConstraintFileID",          StarTagData( GetStarTagCategories().e_DistConstraintList,  "Constraint_file_ID",             util::CPPDataTypes::e_SizeT))),
      e_DistListBlockID(                   AddEnum( "DistListBlockID",                   StarTagData( GetStarTagCategories().e_DistConstraintList,  "Block_ID",                       util::CPPDataTypes::e_SizeT))),
      e_DistListDetails(                   AddEnum( "DistListDetails",                   StarTagData( GetStarTagCategories().e_DistConstraintList,  "Details",                        util::CPPDataTypes::e_String))),

      e_DistValueConstraintID(             AddEnum( "DistValueConstraintID",             StarTagData( GetStarTagCategories().e_DistConstraintValue, "Constraint_ID",                  util::CPPDataTypes::e_SizeT))),
      e_DistValueTreeNodeID(               AddEnum( "DistValueTreeNodeID",               StarTagData( GetStarTagCategories().e_DistConstraintValue, "Tree_node_ID",                   util::CPPDataTypes::e_SizeT))),
      e_DistValueSourceExperimentID(       AddEnum( "DistValueSourceExperimentID",       StarTagData( GetStarTagCategories().e_DistConstraintValue, "Source_experiment_ID",           util::CPPDataTypes::e_SizeT))),
      e_DistValueSpectralPeakID(           AddEnum( "DistValueSpectralPeakID",           StarTagData( GetStarTagCategories().e_DistConstraintValue, "Spectral_peak_ID",               util::CPPDataTypes::e_SizeT))),
      e_DistValueIntensityVal(             AddEnum( "DistValueIntensityVal",             StarTagData( GetStarTagCategories().e_DistConstraintValue, "Intensity_val",                  util::CPPDataTypes::e_Double))),
      e_DistValueIntensityLowerValErr(     AddEnum( "DistValueIntensityLowerValErr",     StarTagData( GetStarTagCategories().e_DistConstraintValue, "Intensity_lower_val_err",        util::CPPDataTypes::e_Double))),
      e_DistValueIntensityUpperValErr(     AddEnum( "DistValueIntensityUpperValErr",     StarTagData( GetStarTagCategories().e_DistConstraintValue, "Intensity_upper_val_err",        util::CPPDataTypes::e_Double))),
      e_DistValueDistanceVal(              AddEnum( "DistValueDistanceVal",              StarTagData( GetStarTagCategories().e_DistConstraintValue, "Distance_val",                   util::CPPDataTypes::e_Double))),
      e_DistValueDistanceLowerBoundVal(    AddEnum( "DistValueDistanceLowerBoundVal",    StarTagData( GetStarTagCategories().e_DistConstraintValue, "Distance_lower_bound_val",       util::CPPDataTypes::e_Double))),
      e_DistValueDistanceUpperBoundVal(    AddEnum( "DistValueDistanceUpperBoundVal",    StarTagData( GetStarTagCategories().e_DistConstraintValue, "Distance_upper_bound_val",       util::CPPDataTypes::e_Double))),
      e_DistValueEntryID(                  AddEnum( "DistValueEntryID",                  StarTagData( GetStarTagCategories().e_DistConstraintValue, "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_DistValueDistanceConstraintListID( AddEnum( "DistValueDistanceConstraintListID", StarTagData( GetStarTagCategories().e_DistConstraintValue, "Distance_constraint_list_ID",    util::CPPDataTypes::e_SizeT))),

      e_RDCID(                             AddEnum( "RDCID",                             StarTagData( GetStarTagCategories().e_RDCConstraint,       "ID",                             util::CPPDataTypes::e_SizeT))),
      e_RDCAssemblyAtomID1(                AddEnum( "RDCAssemblyAtomID1",                StarTagData( GetStarTagCategories().e_RDCConstraint,       "Assembly_atom_ID_1",             util::CPPDataTypes::e_SizeT))),
      e_RDCEntityAssemblyID1(              AddEnum( "RDCEntityAssemblyID1",              StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entity_assembly_ID_1",           util::CPPDataTypes::e_SizeT))),
      e_RDCEntityID1(                      AddEnum( "RDCEntityID1",                      StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entity_ID_1",                    util::CPPDataTypes::e_SizeT))),
      e_RDCCompIndexID1(                   AddEnum( "RDCCompIndexID1",                   StarTagData( GetStarTagCategories().e_RDCConstraint,       "Comp_index_ID_1",                util::CPPDataTypes::e_SizeT))),
      e_RDCSeqID1(                         AddEnum( "RDCSeqID1",                         StarTagData( GetStarTagCategories().e_RDCConstraint,       "Seq_ID_1",                       util::CPPDataTypes::e_SizeT))),
      e_RDCCompID1(                        AddEnum( "RDCCompID1",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Comp_ID_1",                      util::CPPDataTypes::e_String))),
      e_RDCAtomID1(                        AddEnum( "RDCAtomID1",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Atom_ID_1",                      util::CPPDataTypes::e_String))),
      e_RDCResonanceID1(                   AddEnum( "RDCResonanceID1",                   StarTagData( GetStarTagCategories().e_RDCConstraint,       "Resonance_ID_1",                 util::CPPDataTypes::e_SizeT))),
      e_RDCAssemblyAtomID2(                AddEnum( "RDCAssemblyAtomID2",                StarTagData( GetStarTagCategories().e_RDCConstraint,       "Assembly_atom_ID_2",             util::CPPDataTypes::e_SizeT))),
      e_RDCEntityAssemblyID2(              AddEnum( "RDCEntityAssemblyID2",              StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entity_assembly_ID_2",           util::CPPDataTypes::e_SizeT))),
      e_RDCEntityID2(                      AddEnum( "RDCEntityID2",                      StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entity_ID_2",                    util::CPPDataTypes::e_SizeT))),
      e_RDCCompIndexID2(                   AddEnum( "RDCCompIndexID2",                   StarTagData( GetStarTagCategories().e_RDCConstraint,       "Comp_index_ID_2",                util::CPPDataTypes::e_SizeT))),
      e_RDCSeqID2(                         AddEnum( "RDCSeqID2",                         StarTagData( GetStarTagCategories().e_RDCConstraint,       "Seq_ID_2",                       util::CPPDataTypes::e_SizeT))),
      e_RDCCompID2(                        AddEnum( "RDCCompID2",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Comp_ID_2",                      util::CPPDataTypes::e_String))),
      e_RDCAtomID2(                        AddEnum( "RDCAtomID2",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Atom_ID_2",                      util::CPPDataTypes::e_String))),
      e_RDCResonanceID2(                   AddEnum( "RDCResonanceID2",                   StarTagData( GetStarTagCategories().e_RDCConstraint,       "Resonance_ID_2",                 util::CPPDataTypes::e_SizeT))),
      e_RDCVal(                            AddEnum( "RDCVal",                            StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_val",                        util::CPPDataTypes::e_Double))),
      e_RDCLowerBound(                     AddEnum( "RDCLowerBound",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_lower_bound",                util::CPPDataTypes::e_Double))),
      e_RDCUpperBound(                     AddEnum( "RDCUpperBound",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_upper_bound",                util::CPPDataTypes::e_Double))),
      e_RDCValErr(                         AddEnum( "RDCValErr",                         StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_val_err",                    util::CPPDataTypes::e_Double))),
      e_RDCBondLength(                     AddEnum( "RDCBondLength",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_bond_length",                util::CPPDataTypes::e_Double))),
      e_RDCSourceExperimentID(             AddEnum( "RDCSourceExperimentID",             StarTagData( GetStarTagCategories().e_RDCConstraint,       "Source_experiment_ID",           util::CPPDataTypes::e_SizeT))),
      e_RDCAuthAsymID1(                    AddEnum( "RDCAuthAsymID1",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_asym_ID_1",                 util::CPPDataTypes::e_String))),
      e_RDCAuthSeqID1(                     AddEnum( "RDCAuthSeqID1",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_seq_ID_1",                  util::CPPDataTypes::e_String))),
      e_RDCAuthCompID1(                    AddEnum( "RDCAuthCompID1",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_comp_ID_1",                 util::CPPDataTypes::e_String))),
      e_RDCAuthAtomID1(                    AddEnum( "RDCAuthAtomID1",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_atom_ID_1",                 util::CPPDataTypes::e_String))),
      e_RDCAuthAsymID2(                    AddEnum( "RDCAuthAsymID2",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_asym_ID_2",                 util::CPPDataTypes::e_String))),
      e_RDCAuthSeqID2(                     AddEnum( "RDCAuthSeqID2",                     StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_seq_ID_2",                  util::CPPDataTypes::e_String))),
      e_RDCAuthCompID2(                    AddEnum( "RDCAuthCompID2",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_comp_ID_2",                 util::CPPDataTypes::e_String))),
      e_RDCAuthAtomID2(                    AddEnum( "RDCAuthAtomID2",                    StarTagData( GetStarTagCategories().e_RDCConstraint,       "Auth_atom_ID_2",                 util::CPPDataTypes::e_String))),
      e_RDCEntryID(                        AddEnum( "RDCEntryID",                        StarTagData( GetStarTagCategories().e_RDCConstraint,       "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_RDCRDCConstraintListID(            AddEnum( "RDCRDCConstraintListID",            StarTagData( GetStarTagCategories().e_RDCConstraint,       "RDC_constraint_list_ID",         util::CPPDataTypes::e_SizeT))),

      e_RDCListSfCategory(                 AddEnum( "RDCListSfCategory",                 StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Sf_category",                    util::CPPDataTypes::e_String))),
      e_RDCListEntryID(                    AddEnum( "RDCListEntryID",                    StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Entry_ID",                       util::CPPDataTypes::e_String))),
      e_RDCListID(                         AddEnum( "RDCListID",                         StarTagData( GetStarTagCategories().e_RDCConstraintList,   "ID",                             util::CPPDataTypes::e_SizeT))),
      e_RDCListConstraintFileID(           AddEnum( "RDCListConstraintFileID",           StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Constraint_file_ID",             util::CPPDataTypes::e_SizeT))),
      e_RDCListBlockID(                    AddEnum( "RDCListBlockID",                    StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Block_ID",                       util::CPPDataTypes::e_SizeT))),
      e_RDCListDetais(                     AddEnum( "RDCListDetais",                     StarTagData( GetStarTagCategories().e_RDCConstraintList,   "Details",                        util::CPPDataTypes::e_String)))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StarTags::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a tag from a string
    //! @param CATEGORY this tag belongs to
    //! @param DESCRIPTION string description for the enum
    //! @return a tag from a string
    StarTag StarTags::GetTagFromString( const StarTagCategory &CATEGORY, const std::string &DESCRIPTION)
    {
      // iterate through the tag categories
      for
      (
        StarTags::const_iterator tag_itr( GetStarTags().Begin()), tag_itr_end( GetStarTags().End());
        tag_itr != tag_itr_end;
        ++tag_itr
      )
      {
        // if the category has the same category and description
        if( ( *tag_itr)->GetTagCategory() == CATEGORY && ( *tag_itr)->GetDescription() == DESCRIPTION)
        {
          // return it
          return *tag_itr;
        }
      }

      // no matching description found
      return GetStarTags().e_Undefined;
    }

    //! @brief read Star file from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    storage::Map
    <
      StarTagCategory,
      storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
    > StarTags::ReadStarFile( std::istream &ISTREAM)
    {
      storage::Map
      <
        StarTagCategory,
        storage::Pair< storage::Vector< StarTag>, storage::List< std::string> >
      > line_data;

      // read through the file
      while( !ISTREAM.eof())
      {
        // store the line
        std::string buffer;
        std::getline( ISTREAM, buffer);

        // check if it is a start of a data loop
        if( util::TrimString( buffer) == "loop_")
        {
          // get the next line
          std::getline( ISTREAM, buffer);

          // check if the line contains a tag category-tag pair of the format "_Category.Tag"
          storage::Vector< std::string> split_line( util::SplitString( buffer, "."));
          StarTagCategory category
          (
            StarTagCategories::GetCategoryFromString( util::TrimString( split_line.FirstElement()))
          );

          // true if the category is defined
          if( category != GetStarTagCategories().e_Undefined)
          {
            // read in the tag associated with the valid category
            StarTag tag( GetTagFromString( category, util::TrimString( split_line.LastElement())));

            // add the first tag to the vector in the map
            line_data[ category].First().PushBack( tag);

            // read the rest of the tags
            while( !ISTREAM.eof())
            {
              // read in the next line
              std::getline( ISTREAM, buffer);
              split_line = util::SplitString( buffer, ".");

              // if the split line does not have 2 elements, there are no more tags to read in so break
              if( split_line.GetSize() != 2)
              {
                break;
              }

              // add the tag
              line_data[ category].First().PushBack
              (
                GetTagFromString( category, util::TrimString( split_line.LastElement()))
              );
            }

            // iterate through any blank lines
            while( !ISTREAM.eof() && util::TrimString( buffer).empty())
            {
              std::getline( ISTREAM, buffer);
            }

            // iterate through the rows of data that correspond to the tags just read in
            while( !ISTREAM.eof() && util::TrimString( buffer) != "stop_")
            {
              // add the line to the list of strings
              line_data[ category].Second().PushBack( util::TrimString( buffer));
              std::getline( ISTREAM, buffer);
            }
          }

        } // end of data loop
      } // eof
      // end
      return line_data;
    }

    //! @brief gets a chain id from a string since the star format allows authors to identify chains however they want
    //! @param CURRENT_ID current string to be converted
    //! @param PREVIOUS_IDS previously seen strings mapped to chain ids
    //! @return chain id
    char StarTags::GetChainID
    (
      const std::string &CURRENT_ID,
      storage::Map< std::string, char> &PREVIOUS_IDS
    )
    {
      // check the map for the id
      const storage::Map< std::string, char>::const_iterator find_itr( PREVIOUS_IDS.Find( CURRENT_ID));

      // if it is found
      if( find_itr != PREVIOUS_IDS.End())
      {
        // return the corresponding chain
        return find_itr->second;
      }

      // it was not found, so check if the string is 1 character
      if( CURRENT_ID.size() == 1 && CURRENT_ID[ 0] >= 'A' && CURRENT_ID[ 0] <= 'Z')
      {
        // add to map and return
        const char chain_id( CURRENT_ID[ 0]);
        PREVIOUS_IDS[ CURRENT_ID] = chain_id;
        return chain_id;
      }

      // a new, 2+ character id has been passed
      // check if any ids have been set
      if( PREVIOUS_IDS.IsEmpty())
      {
        // set this id to A
        const char initial_id( 'A');
        PREVIOUS_IDS[ CURRENT_ID] = initial_id;
        return initial_id;
      }

      // create a set to find the largets current chain id char
      storage::Set< char> previous_chain_ids;

      // iterate over the map
      for
      (
        storage::Map< std::string, char>::const_iterator map_itr( PREVIOUS_IDS.Begin()),
          map_itr_end( PREVIOUS_IDS.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        previous_chain_ids.Insert( map_itr->second);
      }

      // set the id and return
      char max_chain_id( *previous_chain_ids.ReverseBegin());
      ++max_chain_id;
      PREVIOUS_IDS[ CURRENT_ID] = max_chain_id;
      return max_chain_id;
    }

    //! @brief get access to all star tags
    const StarTags &GetStarTags()
    {
      return StarTags::GetEnums();
    }

  } // namespace nmr

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< nmr::StarTagData, nmr::StarTags>;

  } // namespace util
} // namespace bcl
