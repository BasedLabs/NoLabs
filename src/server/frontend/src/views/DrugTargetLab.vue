<script>
import store from '../storage';
import { api } from '../storage';
import ExperimentsView from './ExperimentsView.vue';
import DrugTargetLabFormView from './DrugTargetLabFormView.vue';
import DrugTargetLabInferenceView from './DrugTargetLabInferenceView.vue';
import DrugTargetList from '../components/AminoAcidLab/DrugTargetList.vue';

export default {
    data() {
        return {
            experimentsState: store.state.drugTarget,
            experimentsApi: api.drugTarget
        }
    },
    components: {
        ExperimentsView,
        DrugTargetLabFormView,
        DrugTargetLabInferenceView,
        DrugTargetList
    }
}
</script>

<template>
    <ExperimentsView :state="experimentsState" :api="experimentsApi">
        <template v-slot:labTitle>
            <h4>Drug target lab</h4>
        </template>
        <template v-slot:labForm="labForm">
            <DrugTargetLabFormView :onFormSubmit="labForm.onFormSubmit"/>
        </template>
        <template v-slot:experimentMetaData="experimentMetaData" :api="experimentsApi">
            <DrugTargetList :experiment="experimentMetaData.experiment" :api="experimentsApi"/>
        </template>
        <template v-slot:lab="lab">
            <DrugTargetLabInferenceView :experiment="lab.experiment"/>
        </template>
    </ExperimentsView >
</template>