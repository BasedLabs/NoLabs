<script>
import store from '../storage';
import { api } from '../storage';
import ExperimentsView from './ExperimentsView.vue';
import ProteinDesignLabFormView from './ProteinDesignLabFormView.vue';
import ProteinDesignLabInferenceView from './ProteinDesignLabInferenceView.vue';

export default {
    data() {
        return {
            experimentsState: store.state.proteinDesign,
            experimentsApi: api.proteinDesign
        }
    },
    components: {
        ExperimentsView,
        ProteinDesignLabFormView,
        ProteinDesignLabInferenceView
    }
}
</script>

<template>
    <ExperimentsView :state="experimentsState" :api="experimentsApi">
        <template v-slot:labTitle>
            <h2>Protein design lab</h2>
        </template>
        <template v-slot:labForm="labForm">
            <ProteinDesignLabFormView :onFormSubmit="labForm.onFormSubmit"/>
        </template>
        <template v-slot:lab="lab">
            <ProteinDesignLabInferenceView :experiment="lab.experiment"/>
        </template>
    </ExperimentsView >
</template>