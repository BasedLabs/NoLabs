<script>
import TargetsList from '../components/AminoAcidLab/TargetsList.vue';
import LigandsList from '../components/AminoAcidLab/LigandList.vue'
import RunDrugTarget from '../components/AminoAcidLab/RunDrugTarget.vue'

export default {
  components: {
    TargetsList,
    LigandsList,
    RunDrugTarget
  },
  props: ['experiment', 'api'],
  data() {
    return {
      currentStep: 1,
      steps: [
        { id: 1, description: 'Select Targets', component: 'TargetsList' },
        { id: 2, description: 'Select Ligands', component: 'LigandsList' },
        { id: 3, description: 'Run Experiment', component: 'RunDrugTarget' }
      ]
    };
  },
  methods: {
    nextStep() {
      if (this.currentStep < this.steps.length) {
        this.currentStep++;
      }
    },
    previousStep() {
      if (this.currentStep > 1) {
        this.currentStep--;
      }
    }
  }
};
</script>

<template>
    <div class="container-fluid">
      <div class="row justify-content-center">
        <div class="col">
          <div class="step-navigation">
            <div class="container-fluid" v-for="step in steps" :key="step.id">
              <div class="row align-items-center mb-3" style="margin-top: 30px;">
                <div class="d-flex align-items-center mb-3">
                    <span v-if="currentStep === step.id" class="step-indicator me-2">&#x1F449;</span>
                    <h4 :class="{ 'text-muted': currentStep !== step.id }" class="mb-0">Step {{ step.id }}: {{ step.description }}</h4>
                </div>
              </div>
              <transition name="slide">
                <component :is="step.component"
                           :experiment="experiment"
                           :api="api"
                           v-if="currentStep === step.id">
                </component>
              </transition>
            </div>
          </div>
          <div class="navigation-buttons my-4">
            <button class="btn btn-secondary mr-2" v-if="currentStep > 1" @click="previousStep">Previous Step</button>
            <button class="btn btn-primary" v-if="currentStep < steps.length" @click="nextStep">Next Step</button>
          </div>
        </div>
      </div>
    </div>
</template>
  
<style>
  .step-navigation {
    display: flex;
    flex-direction: column;
    align-items: flex-start;
    overflow: hidden;
  }
  
  .step-indicator {
    font-size: 1.5rem;
  }
  
  .text-muted {
    color: #6c757d;
  }
  
  .slide-enter-active, .slide-leave-active {
    transition: opacity 0.5s, transform 0.5s;
  }
  
  .slide-enter, .slide-leave-to {
    opacity: 0;
    transform: translateY(20px);
  }
</style>
  