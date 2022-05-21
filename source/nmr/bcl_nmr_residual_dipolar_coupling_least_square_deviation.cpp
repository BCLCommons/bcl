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
