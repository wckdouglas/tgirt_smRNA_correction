import h2o
from h2o.estimators.random_forest import H2ORandomForestEstimator
h2o.init()

class h2o_rf():
    def __init__(self):
        self.rf = H2ORandomForestEstimator()

    def fit(self, X, y):
        '''
        X: pandas dataframe (n x m)
        y: numpy array (n)
        '''
        X_row, X_col  = X.shape
        assert(row > 0 & col > 0)
        x_colnames = X.columns

        assert(y.ndim == 1 & X_row == len(y))
        X['y'] = y.tolist()
        train_df = h2o.H2OFrame.from_python(X)
        self.rf.train(x_colnames, 'y', training_frame=train_df)


    def coeficients(self):
        '''
        return variable importance
        '''
        return self.rf._model_json['output']['variable_importances']\
            .as_data_frame()


    def predict(X):
        '''
        X: dataframe containing training columns
        return: predicted values (list)
        '''
        X = h2o.H2OFrame.from_python(X) 
        y = self.rf.predict(X)
        return y.as_data_frame()['predict'].tolist()

        
    

