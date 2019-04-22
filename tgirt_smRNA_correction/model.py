import h2o
from h2o.estimators.random_forest import H2ORandomForestEstimator
h2o.init()
h2o.no_progress()

class h2o_rf():
    def __init__(self):
        self.rf = H2ORandomForestEstimator(stopping_metric='RMSE')

    def fit(self, X, y):
        '''
        X: pandas dataframe (n x m)
        y: numpy array (n)
        '''
        X_row, X_col  = X.shape
        assert(X_row > 0 and X_col > 0)
        x_colnames = X.columns.tolist()

        assert(y.ndim == 1 and X_row == len(y))
        X.reset_index(inplace=True)
        X['y'] = y.tolist()
        train_df = h2o.H2OFrame.from_python(X)
        self.rf.train(x_colnames, 'y', training_frame=train_df)


    def coeficients(self):
        '''
        return variable importance
        '''
        return self.rf._model_json['output']['variable_importances']\
            .as_data_frame()


    def predict(self, X):
        '''
        X: dataframe containing training columns
        return: predicted values (list)
        '''
        X = h2o.H2OFrame.from_python(X) 
        y = self.rf.predict(X)
        return y.as_data_frame()['predict'].tolist()


    def save_model(self, model_file):
        model_path = h2o.save_model(model=self.rf,  path = model_file)
        return model_path

    def load_model(self, model_path):
        self.rf = h2o.load_model(model_path)
        


        
    

