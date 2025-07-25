#ifndef FWCore_Framework_stream_EDProducerAdaptorBase_h
#define FWCore_Framework_stream_EDProducerAdaptorBase_h
// -*- C++ -*-
//
// Package:     FWCore/Framework
// Class  :     EDProducerAdaptorBase
//
/**\class edm::stream::EDProducerAdaptorBase EDProducerAdaptorBase.h "FWCore/Framework/interface/stream/EDProducerAdaptorBase.h"

 Description: [one line class summary]

 Usage:
    <usage>

*/
//
// Original Author:  Chris Jones
//         Created:  Fri, 02 Aug 2013 18:09:15 GMT
//

// system include files

// user include files
#include "FWCore/Framework/interface/stream/ProducingModuleAdaptorBase.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Provenance/interface/ModuleDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"
#include "FWCore/Concurrency/interface/WaitingTaskHolder.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/RunIndex.h"
#include "FWCore/Utilities/interface/LuminosityBlockIndex.h"

// forward declarations

namespace edm {

  class ModuleCallingContext;
  class ActivityRegistry;

  namespace maker {
    template <typename T>
    class ModuleHolderT;
  }

  namespace stream {
    class EDProducerBase;
    class EDProducerAdaptorBase : public ProducingModuleAdaptorBase<EDProducerBase> {
    public:
      template <typename T>
      friend class edm::maker::ModuleHolderT;
      template <typename T>
      friend class edm::WorkerT;

      EDProducerAdaptorBase();
      EDProducerAdaptorBase(const EDProducerAdaptorBase&) = delete;                   // stop default
      const EDProducerAdaptorBase& operator=(const EDProducerAdaptorBase&) = delete;  // stop default

      // ---------- const member functions ---------------------

      // ---------- static member functions --------------------

      // ---------- member functions ---------------------------

    protected:
      using ProducingModuleAdaptorBase<EDProducerBase>::commit;

    private:
      bool doEvent(EventTransitionInfo const&, ActivityRegistry*, ModuleCallingContext const*);

      void doAcquire(EventTransitionInfo const&, ActivityRegistry*, ModuleCallingContext const*, WaitingTaskHolder&&);

      //For now this is a placeholder
      /*virtual*/ void preActionBeforeRunEventAsync(WaitingTaskHolder,
                                                    ModuleCallingContext const&,
                                                    Principal const&) const noexcept {}
    };
  }  // namespace stream
}  // namespace edm

#endif
