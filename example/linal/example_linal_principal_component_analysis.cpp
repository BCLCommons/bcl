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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "linal/bcl_linal_principal_component_analysis.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_principal_component_analysis.cpp
  //!
  //! @author loweew
  //! @date Sep 6, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalPrincipalComponentAnalysis :
    public ExampleInterface
  {
  public:

    ExampleLinalPrincipalComponentAnalysis *Clone() const
    {
      return new ExampleLinalPrincipalComponentAnalysis( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // creating data set
      float data1[ 9] = { 1.0, 2.0, 3.0, 2.0, 2.0, 4.0, 3.0, 4.0, 3.0};
      linal::Vector< float> vector( 3, 5.0);
      linal::Matrix< float> sqr_mat( 3, 3, 4.0);
      linal::Matrix< float> sym_mat( 3, 3, data1);
      linal::Matrix< float> eigenvec( 3, 3, 0.0f);
      float sort_data[ 9] = { 1.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 2.0};
      float sort_vec_data[ 9] = { 1.0, 3.0, 2.0, 1.0, 3.0, 2.0, 1.0, 3.0, 2.0};

      linal::Matrix< float> eig_val_to_sort( 3, 3, sort_data);
      linal::Matrix< float> eig_vec_to_sort( 3, 3, sort_vec_data);

      storage::Vector< storage::Pair< float, size_t> > sorted_vals_indeces_to_reduce;
      sorted_vals_indeces_to_reduce.PushBack( storage::Pair< float, size_t>( 1.0f, size_t( 2)));
      sorted_vals_indeces_to_reduce.PushBack( storage::Pair< float, size_t>( 2.0f, size_t( 0)));
      sorted_vals_indeces_to_reduce.PushBack( storage::Pair< float, size_t>( 3.0f, size_t( 1)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor
      linal::PrincipalComponentAnalysis< float> pca;

      // clone function
      util::ShPtr< linal::PrincipalComponentAnalysis< float> > pca_clone( pca.Clone());
      BCL_ExampleCheck( pca_clone.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // check operator which returns the sorted eigenvectors v matrix
      linal::Matrix< float> pca_results( pca( sym_mat, 100));
      BCL_ExampleCheckWithinTolerance( double( pca_results( 0, 1)), double( -0.374359), 0.001);
      BCL_ExampleCheckWithinTolerance( double( pca_results( 0, 0)), double( 0.441225), 0.001);

    ////////////////
    // operations //
    ////////////////

      // check SortEigenValuesAndReduce
      storage::VectorND< 2, linal::Matrix< float> > sorted_eigen_vals_and_vecs_reduce( pca.SortEigenValuesAndReduce( eig_val_to_sort, eig_vec_to_sort, 100));
      BCL_ExampleCheck( sorted_eigen_vals_and_vecs_reduce( 0)( 0, 0), float( 3.0));

      // check SortEigenValues
      storage::VectorND< 2, linal::Matrix< float> > sorted_eigen_vals_and_vecs( pca.SortEigenValues( eig_val_to_sort, eig_vec_to_sort));
      BCL_ExampleCheck( sorted_eigen_vals_and_vecs( 0)( 0, 0), float( 3.0));

      // check ReduceEigenValues
      storage::Vector< storage::Pair< float, size_t> > reduced_eig_vals_vecs( pca.ReduceEigenValues( sorted_vals_indeces_to_reduce, 1.0f));
      BCL_ExampleCheck( reduced_eig_vals_vecs( 0).First(), float( 1.0));

      // check ReduceEigenValuesToN
      storage::Vector< storage::Pair< float, size_t> > reduce_eig_vals_vec_to_n( pca.ReduceEigenValuesToN( sorted_vals_indeces_to_reduce, size_t( 2)));
      BCL_ExampleCheck( reduce_eig_vals_vec_to_n( 0).First(), float( 1.0));

      // check GetSortedEigenVectorsValues
      storage::Pair< linal::Matrix< float>, linal::Vector< float> > get_sorted_eig_vec_vals( pca.GetSortedEigenVectorsValues( sym_mat));
      BCL_ExampleCheckWithinTolerance( double( get_sorted_eig_vec_vals.Second()( 0)), double( 8.28822), 0.001);
      BCL_ExampleCheckWithinTolerance( double( get_sorted_eig_vec_vals.First()( 0, 0)), double( 0.441225), 0.001);

      // check ReduceInputMatrix
      pca.ReduceInputMatrix( sym_mat, sym_mat, get_sorted_eig_vec_vals.First(), get_sorted_eig_vec_vals.Second(), 1.0, 100);
      BCL_ExampleCheckWithinTolerance( double( sym_mat( 0, 1)), double( 0.647797), 0.001);
      BCL_ExampleCheckWithinTolerance( double( sym_mat( 0, 0)), double( 3.65697), 0.001);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalPrincipalComponentAnalysis

  const ExampleClass::EnumType ExampleLinalPrincipalComponentAnalysis::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalPrincipalComponentAnalysis())
  );

} // namespace bcl
