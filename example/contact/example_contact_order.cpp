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
#include "contact/bcl_contact_order.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_order.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactOrder :
    public ExampleInterface
  {
  public:

    ExampleContactOrder *Clone() const
    { return new ExampleContactOrder( *this);}

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
      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      //instantiate proteinmodels of chains
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AACaCb, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      contact::Order contact_order_default;

      // from scheme and normalization type
      contact::Order contact_order_absolute( contact::Order::e_Absolute, "co_abs", false);
      contact::Order contact_order_rel_sse( contact::Order::e_RelativeAAsUsed, "co_rel_sse", false);
      contact::Order contact_order_rel_seq( contact::Order::e_RelativeSequenceLength, "co_rel_seq", false);

    /////////////////
    // data access //
    /////////////////

      // scheme
      BCL_ExampleCheck( contact_order_default.GetScheme(), "co");

    ////////////////
    // operations //
    ////////////////

      // contact order of chains
      const double co_chain_abs(     contact_order_absolute.ContactOrder( *model.GetChains().FirstElement()));
      const double co_chain_rel_sse( contact_order_rel_sse.ContactOrder( *model.GetChains().FirstElement()));
      const double co_chain_rel_seq( contact_order_rel_seq.ContactOrder( *model.GetChains().FirstElement()));

      const double co_chain_abs_exp( 27.7778);
      const double co_chain_rel_sse_exp( 0.617284);
      const double co_chain_rel_seq_exp( 0.365497);

      BCL_ExampleIndirectCheckWithinTolerance( co_chain_abs    , co_chain_abs_exp    , 0.001, "absolute contact order of chain");
      BCL_ExampleIndirectCheckWithinTolerance( co_chain_rel_sse, co_chain_rel_sse_exp, 0.001, "relative sse contact order of chain");
      BCL_ExampleIndirectCheckWithinTolerance( co_chain_rel_seq, co_chain_rel_seq_exp, 0.001, "relative sequence contact order of chain");

      // contact order of sequence
      const double co_seq_abs(     contact_order_absolute.ContactOrder( *model.GetChains().FirstElement()->GetSequence()));
      const double co_seq_rel_sse( contact_order_rel_sse.ContactOrder( *model.GetChains().FirstElement()->GetSequence()));
      const double co_seq_rel_seq( contact_order_rel_seq.ContactOrder( *model.GetChains().FirstElement()->GetSequence()));

      const double co_seq_abs_exp( 28.4867);
      const double co_seq_rel_sse_exp( 0.374825);
      const double co_seq_rel_seq_exp( 0.374825);

      BCL_ExampleIndirectCheckWithinTolerance( co_seq_abs    , co_seq_abs_exp    , 0.001, "absolute contact order of seq");
      BCL_ExampleIndirectCheckWithinTolerance( co_seq_rel_sse, co_seq_rel_sse_exp, 0.001, "relative sse contact order of seq");
      BCL_ExampleIndirectCheckWithinTolerance( co_seq_rel_seq, co_seq_rel_seq_exp, 0.001, "relative sequence contact order of seq");

    ////////////////
    // operators //
    ////////////////

      // contact order of protein model
      const double co_pm_abs(     contact_order_absolute( model));
      const double co_pm_rel_sse( contact_order_rel_sse( model));
      const double co_pm_rel_seq( contact_order_rel_seq( model));

      const double co_pm_abs_exp( 27.7778);
      const double co_pm_rel_sse_exp( 0.617284);
      const double co_pm_rel_seq_exp( 0.365497);

      BCL_ExampleIndirectCheckWithinTolerance( co_pm_abs    , co_pm_abs_exp    , 0.001, "absolute contact order of protein model");
      BCL_ExampleIndirectCheckWithinTolerance( co_pm_rel_sse, co_pm_rel_sse_exp, 0.001, "relative sse contact order of protein model");
      BCL_ExampleIndirectCheckWithinTolerance( co_pm_rel_seq, co_pm_rel_seq_exp, 0.001, "relative sequence contact order of protein model");

    //////////////////////
    // input and output //
    //////////////////////

      // write contact order object to file
      WriteBCLObject( contact_order_rel_seq);

      // read
      contact::Order contact_order_read;
      ReadBCLObject( contact_order_read);

      // compare scheme
      BCL_ExampleIndirectCheck( contact_order_read.GetNormalization(), contact_order_rel_seq.GetNormalization(), "write and read bcl object");

      // compare results
      const double co_chain_read( contact_order_read.ContactOrder( *model.GetChains().FirstElement()));
      const double co_seq_read(   contact_order_read.ContactOrder( *model.GetChains().FirstElement()->GetSequence()));
      const double co_pm_read(    contact_order_read( model));

      BCL_ExampleIndirectCheck( co_chain_read, co_chain_rel_seq, "write and read bcl object");
      BCL_ExampleIndirectCheck( co_seq_read  , co_seq_rel_seq  , "write and read bcl object");
      BCL_ExampleIndirectCheck( co_pm_read   , co_pm_rel_seq   , "write and read bcl object");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleContactOrder

  const ExampleClass::EnumType ExampleContactOrder::s_Instance
  (
    GetExamples().AddEnum( ExampleContactOrder())
  );

} // namespace bcl
