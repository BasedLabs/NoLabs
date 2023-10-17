import 'bootswatch/dist/lux/bootstrap.min.css'; // Added this :boom:

import './assets/main.css';
import store from './storage';
import {LoadingPlugin} from 'vue-loading-overlay';
import 'vue-loading-overlay/dist/css/index.css';

import { createApp } from 'vue';
import App from './App.vue';
import router from './router';

const app = createApp(App);

app.use(LoadingPlugin);
app.use(router);
app.use(store);

app.mount('#app');
